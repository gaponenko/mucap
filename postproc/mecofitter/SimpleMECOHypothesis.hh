// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_SimpleMECOHypothesis_hh
#define postproc_mecofitter_SimpleMECOHypothesis_hh

#include <functional>

#include "MuCapUtilities/inc/mecoSpectrum.hh"

namespace mucap {
  namespace fitter {

    class SimpleMECOHypothesis : public std::unary_function<double,double> {
        MECOPars pars_;

    public:

      double operator()(double ek) const;

      const MECOPars& pars() const { return pars_; }

      explicit SimpleMECOHypothesis(const MECOPars& p);

      // This ugly ctr is for use with the ROOT interface
      explicit SimpleMECOHypothesis(const double *x);
    };


  } // fitter
} // mucap

#endif/*postproc_mecofitter_SimpleMECOHypothesis_hh*/
