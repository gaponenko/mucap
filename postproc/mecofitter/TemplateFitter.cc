#include "TemplateFitter.hh"

#include <iostream>
#include <numeric>

#include "cetlib/exception.h"

#include "getHisto.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "Math/Factory.h"

#include "MuCapUtilities/inc/mecoSpectrum.hh"
#include "SimpleMECOHypothesis.hh"
#include "applyHypothesis.hh"

namespace mucap {
  namespace fitter {

    using std::string;
    using fhicl::ParameterSet;

    //================================================================
    TemplateFitter::TemplateFitter(const fhicl::ParameterSet& pset)
      : pset_(pset)
      , origData_(readDataHisto(pset.get<ParameterSet>("data_inputs")))
      , origMC_(readMCHisto(pset.get<ParameterSet>("mc_inputs")))
      , region_(computeFitRegion(pset_.get<ParameterSet>("cuts"), origData_.get(), origMC_.get()))

      , func_(extractFitData(origData_.get(), region_),
              extractFitMC(origMC_.get(), region_))

      , min_(createMinimizer(pset.get<ParameterSet>("minimizer")))
    {
      MECOPars fg; // first guess parameters

      fg.A = 1.;  // will estimate the normalization below
      fg.Tth = 1.4;
      fg.alpha = 1.33;
      // We fit spectrum with T0=3.1 using MC data with T0=3.5
      // 3.1*3.5/(3.5-3.1) = 27.125
      fg.T0 = 27.125;

      //----------------
      // Determine the initial normalization:

      const FitData& data = func_.getFitData();
      const double dataIntegral = std::accumulate(data.begin(), data.end(), 0.);

      SimpleMECOHypothesis testSpectrum(fg);
      const FitData mc = applyHypothesis(func_.getFitMC(), testSpectrum);
      const double mcIntegral = std::accumulate(mc.begin(), mc.end(), 0.);

      fg.A = dataIntegral/mcIntegral;

      //----------------

      min_->SetFunction(func_);

      min_->SetVariable(0, "A",     fg.A,     0.1);
      min_->SetVariable(1, "Tth",   fg.Tth,   0.1);
      min_->SetVariable(2, "alpha", fg.alpha, 0.1);
      min_->SetVariable(3, "T0",    fg.T0,    0.1*fg.T0);

      min_->SetMaxFunctionCalls(1000000);
      min_->SetMaxIterations(100000);
      min_->SetTolerance(0.001);
    }

    //================================================================
    TH2 *TemplateFitter::readDataHisto(const ParameterSet& data_inputs) {
      TH2 * res = dynamic_cast<TH2*>(getHisto(data_inputs.get<string>("filename"),
                                              data_inputs.get<string>("histogram")));
      if(!res) {
        throw cet::exception("INPUTS")<<"The data histo is not a TH2\n";
      }

      return res;
    }

    //================================================================
    TH3 *TemplateFitter::readMCHisto(const ParameterSet& mc_inputs) {
      TH3 * res = dynamic_cast<TH3*>(getHisto(mc_inputs.get<string>("filename"),
                                              mc_inputs.get<string>("histogram")));
      if(!res) {
        throw cet::exception("INPUTS")<<"The MC histo is not a TH3\n";
      }

      return res;
    }

    //================================================================
    FitRegion TemplateFitter::computeFitRegion(const fhicl::ParameterSet& cuts,
                                               const TH2 *hdata,
                                               const TH3 *hmc)
    {
      const TH2* mcproj = dynamic_cast<TH2*>(hmc->Project3D("yx"));
      if(!mcproj) {
        throw cet::exception("INPUTS")<<"Error when projecting MC input\n";
      }

      FitRegion res;
      const unsigned minBinContent = cuts.get<unsigned>("minProjectedMCBinContent");
      int minPlanes = cuts.get<unsigned>("minPlanes");
      int maxPlanes = cuts.get<unsigned>("maxPlanes");
      int minCells  = cuts.get<unsigned>("minCells");
      int maxCells  = cuts.get<unsigned>("maxCells");

      // We assume 1 cell per X bin, 1 plane per Y bin,
      // starting at (0 cells, 0 planes) == ROOT bin (1,1)
      if(mcproj->GetNbinsX() < maxCells + 1) {
        throw cet::exception("INPUTS")
          <<"maxCells = "<<maxCells<<" exceeds the range of the input MC histogram: "
          <<mcproj->GetNbinsX()
          <<"\n";
      }
      if(mcproj->GetNbinsY() < maxPlanes + 1) {
        throw cet::exception("INPUTS")
          <<"maxPlanes = "<<maxPlanes<<" exceeds the range of the input MC histogram\n";
      }

      if(hdata->GetNbinsX() < maxCells + 1) {
        throw cet::exception("INPUTS")
          <<"maxCells = "<<maxCells<<" exceeds the range of the input data histogram\n";
      }
      if(hdata->GetNbinsY() < maxPlanes + 1) {
        throw cet::exception("INPUTS")
          <<"maxPlanes = "<<maxPlanes<<" exceeds the range of the input data histogram\n";
      }

      for(int ix = 1+minCells; ix <= 1+maxCells; ++ix) {
        for(int iy = 1+minPlanes; iy <= 1+maxPlanes; ++iy) {
          const double nmc = mcproj->GetBinContent(ix, iy);
          if(minBinContent <= nmc) {
            res.push_back(Bin(ix,iy));
          }
          else { // We are about to exclude the bin based on low MC statistics
            // Check that it has no data entries
            if(hdata->GetBinContent(ix,iy) > 0.) {
              // throw cet::exception("INPUTS")
              std::cout<<"WARNING: "
                <<"Bin("<<ix<<", "<<iy<<") "
                <<"has "<<hdata->GetBinContent(ix,iy)
                <<" data entries but insufficient MC statistics = "
                <<mcproj->GetBinContent(ix, iy)<<"\n";
            }
          }
        }
      }

      return res;
    }

    //================================================================
    FitData TemplateFitter::extractFitData(const TH2* hh, const FitRegion& r) {
      FitData res;

      res.reserve(r.size());
      for(FitRegion::const_iterator i=r.begin(); i!=r.end(); ++i) {
        res.push_back(hh->GetBinContent(i->ix, i->iy));
      }

      return res;
    }

    //================================================================
    FitMC TemplateFitter::extractFitMC(const TH3* hh, const FitRegion& r) {
      FitMC res;

      const TAxis *az = hh->GetZaxis();
      const int nzbins = az->GetNbins();
      res.reserve(nzbins);
      for(int iz=1; iz<=nzbins; ++iz) {
        FitMCSlice slice;
        slice.energy = az->GetBinCenter(iz);
        slice.values.reserve(r.size());
        for(FitRegion::const_iterator i=r.begin(); i!=r.end(); ++i) {
          slice.values.push_back(hh->GetBinContent(i->ix, i->iy, iz));
        }
        res.push_back(slice);
      }

      return res;
    }

    //================================================================
    ROOT::Math::Minimizer* TemplateFitter::createMinimizer(const ParameterSet &pset) {
      return ROOT::Math::Factory::
        CreateMinimizer(pset.get<string>("type"), pset.get<string>("strategy"));
    }

    //================================================================
    TH2* TemplateFitter::visualizeFitRegion(const TH2 *binningTemplate,
                                            const FitRegion &r) const
    {
      const TAxis *ax = binningTemplate->GetXaxis();
      const TAxis *ay = binningTemplate->GetYaxis();

      TH2D *res = new TH2D("fitRegion", "fitRegion",
                           ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                           ay->GetNbins(), ay->GetXmin(), ay->GetXmax()
                           );

      for(FitRegion::const_iterator i=r.begin(); i!=r.end(); ++i) {
        res->SetBinContent(i->ix, i->iy, 1.);
      }

      return res;
    }

    //================================================================
    TH2* TemplateFitter::visualizeFitData(const std::string& hname,
                                          const FitData& data,
                                          const TH2 *binningTemplate,
                                          const FitRegion &r
                                          ) const
    {
      const TAxis *ax = binningTemplate->GetXaxis();
      const TAxis *ay = binningTemplate->GetYaxis();

      TH2D *res = new TH2D(hname.c_str(), hname.c_str(),
                           ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                           ay->GetNbins(), ay->GetXmin(), ay->GetXmax()
                           );

      for(unsigned i=0; i<r.size(); ++i) {
        res->SetBinContent(r[i].ix, r[i].iy, data[i]);
      }

      return res;
    }

    //================================================================
  }
}
