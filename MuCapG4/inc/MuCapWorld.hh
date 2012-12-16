// Andrei Gaponenko, 2012, based on Mu2eStudyWorld by Krzysztof Genser.

#ifndef MuCapG4_inc_MuCapWorld_hh
#define MuCapG4_inc_MuCapWorld_hh

#include "fhiclcpp/ParameterSet.h"

#include "Mu2eG4/inc/Mu2eUniverse.hh"
class G4VPhysicalVolume;

namespace mu2e {
  class VolumeInfo;

  class MuCapWorld : public Mu2eUniverse {
    fhicl::ParameterSet geom_;

    bool forceAuxEdgeVisible_;
    bool doSurfaceCheck_;
    bool placePV_;

    static const std::string targetModuleType_;

    // Returns the number of planes in the module
    unsigned constructChamberModule(unsigned moduleNumber,
                                    unsigned planeNumberOffset,
                                    const fhicl::ParameterSet& modulePars,
                                    const VolumeInfo& parent);

    void constructFoil(unsigned imodule,
                       unsigned ifoil,
                       const fhicl::ParameterSet& foilPars,
                       const CLHEP::Hep3Vector& centerInParent,
                       const VolumeInfo& parent);

  public:

    explicit MuCapWorld(const fhicl::ParameterSet& pset);
    ~MuCapWorld();

    // Construct everything.
    // The non-const return type is eventually required
    // by G4VUserDetectorConstruction::Construct();
    G4VPhysicalVolume * construct();
  };

} // end namespace mu2e

#endif /* MuCapG4_inc_MuCapWorld_hh */
