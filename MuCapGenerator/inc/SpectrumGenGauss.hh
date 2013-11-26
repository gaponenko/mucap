// An exponentially falling spectrum ~ exp(-x/scale) in a [min, max] range.
//
// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_SpectrumGenGauss_hh
#define MuCapGenerator_inc_SpectrumGenGauss_hh

#include "CLHEP/Random/RandGaussQ.h"

#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"

namespace mucap {

  class SpectrumGenGauss: public ISpectrumGenerator {
  public:

    virtual double generate();

    SpectrumGenGauss(double mean, double sigma, CLHEP::HepRandomEngine& rng);

  private:
    CLHEP::RandGaussQ gauss_;
  };
}

#endif/*MuCapGenerator_inc_SpectrumGenGauss_hh*/
