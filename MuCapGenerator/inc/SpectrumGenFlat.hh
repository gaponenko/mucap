// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_SpectrumGenFlat_hh
#define MuCapGenerator_inc_SpectrumGenFlat_hh

#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"
#include "CLHEP/Random/RandFlat.h"

namespace mucap {

  class SpectrumGenFlat: public ISpectrumGenerator {
  public:

    virtual double generate();

    SpectrumGenFlat(double center,
                    double halfWidth,
                    CLHEP::HepRandomEngine& rng);

  private:
    CLHEP::RandFlat flat_;
    double center_;
    double halfWidth_;
  };
}

#endif/*MuCapGenerator_inc_SpectrumGenFlat_hh*/
