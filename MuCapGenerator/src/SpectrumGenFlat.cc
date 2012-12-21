// Andrei Gaponenko, 2012

#include "MuCapGenerator/inc/SpectrumGenFlat.hh"

#include <cmath>

namespace mucap {

  SpectrumGenFlat::SpectrumGenFlat(double c, double hw, CLHEP::HepRandomEngine& rng)
    : flat_(rng)
    , center_(c)
    , halfWidth_(hw)
  {}

  double SpectrumGenFlat::generate() {
    return center_ + 2*halfWidth_*(flat_.fire() - 0.5);
  }

}
