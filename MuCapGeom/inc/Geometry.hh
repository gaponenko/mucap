// Andrei Gaponenko, 2012

#ifndef MuCapGeom_inc_Geometry_hh
#define MuCapGeom_inc_Geometry_hh

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {

  enum class WPType { PC, DC };
  struct WPInfo { WPType wpt; unsigned int localPlane; };

  class Geometry {
  public:

    Geometry(const fhicl::ParameterSet&, art::ActivityRegistry&);

    const fhicl::ParameterSet& pset() const { return pset_; }

    const CLHEP::Hep3Vector& targetCenter() const { return targetCenter_; }
    double targetThickness() const { return targetThickness_; }

    // globalPlaneNumber is 1-based
    const WPInfo& wpByGlobalNumber(unsigned globalPlaneNumber) const { return wpi_.at(globalPlaneNumber-1); }

  private:
    fhicl::ParameterSet pset_;

    CLHEP::Hep3Vector targetCenter_;
    double targetThickness_;

    std::vector<WPInfo> wpi_;
  };
}

DECLARE_ART_SERVICE(mucap::Geometry, LEGACY)
#endif/*MuCapGeom_inc_Geometry_hh*/
