// A proton spectrum shape used in MECO note 034 (Mu2e docdb-105).
//
// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_SpectrumGenMECO_hh
#define MuCapGenerator_inc_SpectrumGenMECO_hh

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"
#include "MuCapUtilities/inc/mecoSpectrum.hh"

namespace mucap {

  class SpectrumGenMECO: public ISpectrumGenerator {
  public:

    virtual double generate();

    SpectrumGenMECO(const MECOPars& pars,
                    CLHEP::HepRandomEngine& rng,
                    unsigned maxIter);

  private:
    CLHEP::RandExponential exp_;
    CLHEP::RandFlat flat_;
    MECOPars pars_;
    unsigned maxIter_;
  };
}

#endif/*MuCapGenerator_inc_SpectrumGenMECO_hh*/
