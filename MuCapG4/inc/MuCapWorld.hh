// Andrei Gaponenko, 2012, based on Mu2eStudyWorld by Krzysztof Genser.

#ifndef MuCapG4_inc_MuCapWorld_hh
#define MuCapG4_inc_MuCapWorld_hh

#include "MuCapGeom/inc/Geometry.hh"

#include "fhiclcpp/ParameterSet.h"

#include "Mu2eG4/inc/Mu2eUniverse.hh"

class G4VPhysicalVolume;

namespace mu2e { class VolumeInfo; }

namespace mucap {

  class MuCapWorld : public mu2e::Mu2eUniverse {
    const Geometry *geom_; // non-owning pointer

    bool forceAuxEdgeVisible_;
    bool doSurfaceCheck_;
    bool placePV_;

    static const std::string targetModuleType_;

    mu2e::VolumeInfo constructDetectorEnclosure(const mu2e::VolumeInfo& parent, const fhicl::ParameterSet& pset);

    // Returns the number of planes in the module
    unsigned constructChamberModule(unsigned moduleNumber,
                                    unsigned planeNumberOffset,
                                    const fhicl::ParameterSet& modulePars,
                                    const mu2e::VolumeInfo& parent);

    void constructFoil(unsigned imodule,
                       unsigned ifoil,
                       const fhicl::ParameterSet& foilPars,
                       const CLHEP::Hep3Vector& centerInParent,
                       const mu2e::VolumeInfo& parent);

    void constructCathodeSupport(unsigned imodule,
                                 unsigned ifoil,
                                 const fhicl::ParameterSet& supportPars,
                                 double rIn, double rOut,
                                 const CLHEP::Hep3Vector& centerInParent,
                                 const mu2e::VolumeInfo& parent);

    void constructGlassFrame(unsigned imodule,
                             unsigned iframe,
                             const fhicl::ParameterSet& framePars,
                             const CLHEP::Hep3Vector& centerInParent,
                             const mu2e::VolumeInfo& parent);


    void constructDriftPlane(unsigned globalPlaneNumber,
                             const fhicl::ParameterSet& detail,
                             double planeRotation, // radians
                             double zwire, // in parent
                             double driftZmin, // in parent
                             double driftZmax, // in parent
                             const mu2e::VolumeInfo& parent
                             );
  public:

    explicit MuCapWorld(const Geometry& geom);

    ~MuCapWorld();

    // Construct everything.
    // The non-const return type is eventually required
    // by G4VUserDetectorConstruction::Construct();
    G4VPhysicalVolume * construct();
  };

} // end namespace mucap

#endif /* MuCapG4_inc_MuCapWorld_hh */
