// Generic cylinder along the z axis.
//
// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_AngleGenUniform_hh
#define MuCapGenerator_inc_AngleGenUniform_hh

#include "MuCapInterfaces/inc/IAngleGenerator.hh"

#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

namespace mucap {


  class AngleGenUniform: public IAngleGenerator {
  public:

    virtual CLHEP::Hep3Vector generate() {
      return sphere_.fire();
    }

    AngleGenUniform(const mu2e::RandomUnitSphereParams& angles, CLHEP::HepRandomEngine& rng)
      : sphere_(rng, angles)
    {}

  private:

    mu2e::RandomUnitSphere sphere_;

  };
}

#endif/*MuCapGenerator_inc_AngleGenUniform_hh*/
