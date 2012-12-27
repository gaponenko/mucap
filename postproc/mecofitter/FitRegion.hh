// FitRegion is a list of bins (ix,iy) of a two-dimensional
// histogram that should be included in the fit.
//
// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_FitRegion_hh
#define postproc_mecofitter_FitRegion_hh

#include <vector>

#include "fhiclcpp/ParameterSet.h"

namespace mucap {
  namespace fitter {

    struct Bin {
      int ix;
      int iy;
      Bin(int a, int b) : ix(a), iy(b) {}
    };

    typedef std::vector<Bin> FitRegion;

  } // fitter
} // mucap

#endif/*postproc_mecofitter_FitRegion_hh*/
