#include "SimpleMECOHypothesis.hh"

namespace mucap {
  namespace fitter {

    double SimpleMECOHypothesis::operator()(double ek) const {
      return mecoSpectrum(ek, pars_);
    }

    SimpleMECOHypothesis::SimpleMECOHypothesis(const MECOPars& p)
      : pars_(p)
    {}

    SimpleMECOHypothesis::SimpleMECOHypothesis(const double *x) {
      pars_.A = x[0];
      pars_.Tth = x[1];
      pars_.alpha = x[2];
      pars_.T0 = x[3];
    }

  }
}
