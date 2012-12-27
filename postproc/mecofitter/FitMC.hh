// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_FitMC_hh
#define postproc_mecofitter_FitMC_hh

#include <vector>

#include "FitData.hh"

namespace mucap {
  namespace fitter {

    struct FitMCSlice {
      double energy;
      FitData values;
      FitMCSlice() : energy(), values() {}
    };

    typedef std::vector<FitMCSlice> FitMC;

  } // fitter
} // mucap

#endif/*postproc_mecofitter_FitMC_hh*/
