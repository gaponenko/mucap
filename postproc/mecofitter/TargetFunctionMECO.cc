#include "TargetFunctionMECO.hh"

#include <cassert>

#include "MuCapUtilities/inc/mecoSpectrum.hh"
#include "SimpleMECOHypothesis.hh"
#include "applyHypothesis.hh"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mucap {
  namespace fitter {

    //================================================================
    TargetFunctionMECO::TargetFunctionMECO(const FitData& data, const FitMC& mc)
      : data_(data), mc_(mc)
    {
      assert(!mc_.empty());
      assert(data_.size() == mc_.front().values.size());
      assert(data_.size());
    }

    //================================================================
    ROOT::Math::IBaseFunctionMultiDim* TargetFunctionMECO::Clone() const {
      return new TargetFunctionMECO(*this);
    }

    //================================================================
    unsigned int TargetFunctionMECO::NDim() const {
      return 4; // the number of variables in MECOPars
    }

    //================================================================
    // Compute -2*ln(lambda), S.Baker, R.D.Cousins, NIMA221(1984)437
    double TargetFunctionMECO::DoEval(const double* x) const {

      SimpleMECOHypothesis sm(x);

      FitData prediction(applyHypothesis(mc_, sm));

      double res(0);
      for(unsigned i=0; i<data_.size(); ++i) {
        if(prediction[i] > 0.) {
          res += prediction[i] - data_[i]
            + ( data_[i] > 0 ?  data_[i] * log(data_[i]/prediction[i]) : 0.);
        }
        else {
          if(data_[i] > 0.) {
            throw cet::exception("RUNTIME")<<__func__
                                           <<": got bin with prediction="<<prediction[i]
                                           <<" but data= "<<data_[i]<<" > 0: i = "<<i
                                           <<" at "<<sm.pars()
                                           <<"\n";
          }
        }
      }

      return 2*res;
    }

    //================================================================
  }
}
