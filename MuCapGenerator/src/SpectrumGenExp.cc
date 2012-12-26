// Andrei Gaponenko, 2012


#include "MuCapGenerator/inc/SpectrumGenExp.hh"

#include <cmath>

#include "fhiclcpp/ParameterSet.h"


#include <iostream>
namespace mucap {

  SpectrumGenExp::SpectrumGenExp(double scale, double min, double max, CLHEP::HepRandomEngine& rng)
    : flat_(rng)
    , scale_(scale)
    , r1_(exp(-max/scale))
    , r2_(exp(-min/scale))
  {}

  double SpectrumGenExp::generate() {
    return -scale_ * log(flat_.fire(r1_, r2_));
  }

} // namespace mucap
