// Andrei Gaponenko, 2012

#include "MuCapGenerator/inc/PosGenCylinder.hh"

#include <cmath>

namespace mucap {

  PosGenCylinder::PosGenCylinder(const CLHEP::Hep3Vector center,
                                 double radius,
                                 double halfdz,
                                 CLHEP::HepRandomEngine& rng)
    : flat_(rng)
    , center_(center)
    , radius_(radius)
    , halfdz_(halfdz)
  {}

  CLHEP::Hep3Vector PosGenCylinder::generate() {
    static const double pi = 4.*atan(1.);

    const double z = 2*halfdz_*(flat_.fire() - 0.5);
    const double r =  radius_ * sqrt(flat_.fire());
    const double phi = 2*pi*flat_.fire();

    return center_ + CLHEP::Hep3Vector(r*cos(phi), r*sin(phi), z);
  }

}
