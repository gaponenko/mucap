// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_TargetFunctionMECO_hh
#define postproc_mecofitter_TargetFunctionMECO_hh

#include "Math/IFunction.h"

#include "FitRegion.hh"
#include "FitData.hh"
#include "FitMC.hh"

namespace mucap {
  namespace fitter {


    class TargetFunctionMECO : public ROOT::Math::IBaseFunctionMultiDim {
      FitData data_;
      FitMC   mc_;

    public:
      // Methods required by the interface
      virtual IBaseFunctionMultiDim * Clone() const;

      virtual unsigned int NDim() const;

      virtual double DoEval(const double* x) const;

      // Extra stuff
      TargetFunctionMECO(const FitData& data, const FitMC& mc);

      const FitData& getFitData() const { return data_; }
      const FitMC&   getFitMC() const   { return mc_; }
    };

  } // fitter
} // mucap

#endif/*postproc_mecofitter_TargetFunctionMECO_hh*/
