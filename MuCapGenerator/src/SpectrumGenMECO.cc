// Andrei Gaponenko, 2012


#include "MuCapGenerator/inc/SpectrumGenMECO.hh"

#include <cmath>

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"

namespace mucap {

  SpectrumGenMECO::SpectrumGenMECO(const MECOPars& pars,
                                   CLHEP::HepRandomEngine& rng,
                                   unsigned maxIter)
    : exp_(rng)
    , flat_(rng)
    , pars_(pars)
    , maxIter_(maxIter)
  {}

  double SpectrumGenMECO::generate() {

    for(unsigned i=0; i<maxIter_; ++i) {

      // The exponential part of the spectrum is generated analytically (via CLHEP)
      // then accept/reject is done on the remaining factor.
      const double ek = pars_.Tth + exp_.fire(1./pars_.T0inv);

      const double factor = std::pow(1 - pars_.Tth/ek, pars_.alpha);
      const double r = flat_.fire();

      if(r < factor) {
        return ek;
      }
    }

    throw cet::exception("RUNTIME")<<__func__<<": too many iteration, the limit is "<<maxIter_<<"\n";
  }

} // namespace mucap
