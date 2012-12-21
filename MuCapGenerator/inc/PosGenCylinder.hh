// Generic cylinder along the z axis.
//
// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_PosGenCylinder_hh
#define MuCapGenerator_inc_PosGenCylinder_hh

#include "MuCapInterfaces/inc/IPositionGenerator.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {


  class PosGenCylinder: public IPositionGenerator {
  public:
    virtual CLHEP::Hep3Vector generate();

    PosGenCylinder(const CLHEP::Hep3Vector center,
                   double radius,
                   double halfdz,
                   CLHEP::HepRandomEngine& rng);

  private:
    CLHEP::RandFlat flat_;
    CLHEP::Hep3Vector center_;
    double radius_;
    double halfdz_;
  };
}

#endif/*MuCapGenerator_inc_PosGenCylinder_hh*/
