// Andrei Gaponenko, 2012

#ifndef MuCapGeom_inc_Geometry_hh
#define MuCapGeom_inc_Geometry_hh

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"

#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {
  class Geometry {
  public:

    Geometry(const fhicl::ParameterSet&, art::ActivityRegistry&);

    const fhicl::ParameterSet& pset() const { return pset_; }

    const CLHEP::Hep3Vector& targetCenter() const { return targetCenter_; }
    double targetThickness() const { return targetThickness_; }

  private:
    fhicl::ParameterSet pset_;

    CLHEP::Hep3Vector targetCenter_;
    double targetThickness_;
  };
}

#endif/*MuCapGeom_inc_Geometry_hh*/
