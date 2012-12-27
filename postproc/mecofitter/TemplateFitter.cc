#include "TemplateFitter.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "getHisto.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mucap {
  namespace fitter {

    using std::string;
    using fhicl::ParameterSet;

    //================================================================
    TemplateFitter::TemplateFitter(const fhicl::ParameterSet& pset)
      : pset_(pset)
    {
      const ParameterSet mc_inputs(pset_.get<ParameterSet>("mc_inputs"));
      TH1* h1 = getHisto(mc_inputs.get<string>("filename"), mc_inputs.get<string>("histogram"));
      AGDEBUG("hh = "<<h1);

      const TH3* hmc = dynamic_cast<TH3*>(h1);
      if(!hmc) {
        throw cet::exception("INPUTS")
          <<"MC input histogram is not a TH3\n";
      }

      const TH2* hmc2 = dynamic_cast<TH2*>(hmc->Project3D("yx"));
      if(!hmc2) {
        throw cet::exception("INPUTS")
          <<"Error when projecting MC input\n";
      }

      const ParameterSet data_inputs(pset_.get<ParameterSet>("data_inputs"));
      const TH2* hdata = dynamic_cast<TH2*>(getHisto(data_inputs.get<string>("filename"),
                                                     data_inputs.get<string>("histogram")));

      region_ = computeFitRegion(pset_.get<ParameterSet>("cuts"), hdata, hmc2);

      FitData dd = getFitData(hdata, region_);
      FitMC   mm = getFitMC(hmc, region_);

      //----------------------------------------------------------------
      const std::string outfilename(pset.get<ParameterSet>("outputs").get<string>("filename"));
      TFile out(outfilename.c_str(), "RECREATE");
      if(!out.IsOpen()) {
        throw cet::exception("INPUTS")
          <<"Error opening output file "<<outfilename<<"\n";
      }

      visualizeFitRegion(hmc2, region_)->Write();

      hmc2->Write();

      out.Close();
    }

    //================================================================
    FitRegion TemplateFitter::computeFitRegion(const fhicl::ParameterSet& cuts,
                                               const TH2 *hdata,
                                               const TH2 *mcproj) const
    {
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
    FitData TemplateFitter::getFitData(const TH2* hh, const FitRegion& r) const {
      FitData res;

      res.reserve(r.size());
      for(FitRegion::const_iterator i=r.begin(); i!=r.end(); ++i) {
        res.push_back(hh->GetBinContent(i->ix, i->iy));
      }

      return res;
    }

    //================================================================
    FitMC TemplateFitter::getFitMC(const TH3* hh, const FitRegion& r) const {
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
  }
}
