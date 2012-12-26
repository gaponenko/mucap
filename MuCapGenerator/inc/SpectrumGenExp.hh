// An exponentially falling spectrum ~ exp(-x/scale) in a [min, max] range.
//
// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_SpectrumGenExp_hh
#define MuCapGenerator_inc_SpectrumGenExp_hh

#include "CLHEP/Random/RandFlat.h"

#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"

namespace mucap {

  class SpectrumGenExp: public ISpectrumGenerator {
  public:

    virtual double generate();

    SpectrumGenExp(double scale, double min, double max, CLHEP::HepRandomEngine& rng);

  private:
    CLHEP::RandFlat flat_;
    double scale_;
    double r1_;
    double r2_;
  };
}

#endif/*MuCapGenerator_inc_SpectrumGenExp_hh*/
