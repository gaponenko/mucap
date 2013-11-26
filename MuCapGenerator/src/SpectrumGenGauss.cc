// Andrei Gaponenko, 2012


#include "MuCapGenerator/inc/SpectrumGenGauss.hh"

#include <cmath>

#include "fhiclcpp/ParameterSet.h"


#include <iostream>
namespace mucap {

  SpectrumGenGauss::SpectrumGenGauss(double mean, double sigma, CLHEP::HepRandomEngine& rng)
    : gauss_(rng, mean, sigma)
  {}

  double SpectrumGenGauss::generate() {
    return gauss_.fire();
  }

} // namespace mucap
