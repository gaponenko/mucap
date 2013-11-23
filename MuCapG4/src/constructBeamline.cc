// Andrei Gaponenko, 2013

#include "MuCapG4/inc/MuCapWorld.hh"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <utility>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4ClassicalRK4.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4GDMLParser.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"

#include "MuCapG4/inc/MuCapSD.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mucap {

  using namespace std;
  using namespace mu2e;
  using fhicl::ParameterSet;

  void MuCapWorld::constructBeamline(const VolumeInfo& parent, const fhicl::ParameterSet& pset) {
    const double zmin = pset.get<double>("zmin");
    const double zend = pset.get<double>("gabs_zend");
    const double gabs_length = pset.get<double>("gabs_outer_length");
    const double gabs_rOut = pset.get<double>("gabs_outer_radius");

    TubsParams vac_params(0., gabs_rOut, (zend - zmin)/2);
    const CLHEP::Hep3Vector vacCenterInWorld(0,0, (zend + zmin)/2);

    const VolumeInfo vac = nestTubs("BEAMVAC",
                                 vac_params,
                                 findMaterialOrThrow(pset.get<string>("vacuum_material")),
                                 0, // no rotation
                                 vacCenterInWorld - parent.centerInWorld,
                                 parent,
                                 0,
                                 pset.get<bool>("vacuum_visible"),
                                 G4Colour::White(),
                                 pset.get<bool>("vacuum_solid"),
                                 forceAuxEdgeVisible_,
                                 placePV_,
                                 doSurfaceCheck_
                                 );

    // create the beam pipe wall that extends to the degrader
    TubsParams bpipe_params(gabs_rOut - pset.get<double>("wall_thickness"), gabs_rOut, (zend - gabs_length - zmin)/2);
    const CLHEP::Hep3Vector bpipeCenterInWorld(0,0, (zend - gabs_length + zmin)/2);
    nestTubs("BeamPipeWall",
             bpipe_params,
             findMaterialOrThrow(pset.get<string>("wall_material")),
             0, // no rotation
             bpipeCenterInWorld - vac.centerInWorld,
             vac,
             0,
             pset.get<bool>("wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );
    
    constructGasDegrader(vac, zend, gabs_rOut, gabs_length, pset.get<ParameterSet>("gabs"));
  }

  //================================================================
  void MuCapWorld::constructGasDegrader(const VolumeInfo& parent, double zend, double rOut, double lOut, const fhicl::ParameterSet& pset) {
    const double rIn = rOut - pset.get<double>("wall_thickness");
    const double exitFoilThick = pset.get<double>("exit_foil_thickness");
    const double foilDistance = pset.get<double>("foilToFoilNominalDistance");
    const double rVacWin = pset.get<double>("vacuum_window_radius");
    const double tVacWin = pset.get<double>("vacuum_window_thickness");
    const double vacHolderThick = pset.get<double>("vacuum_holder_thickness");

    nestTubs("GabsGasWall",
             TubsParams(rIn, rOut, lOut/2),
             findMaterialOrThrow(pset.get<string>("wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, zend - lOut/2) - parent.centerInWorld,
             parent,
             0,
             pset.get<bool>("wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );
    
    nestTubs("GabsExitFoil",
             TubsParams(0., rIn, exitFoilThick/2),
             findMaterialOrThrow(pset.get<string>("exit_foil_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, zend - exitFoilThick/2) - parent.centerInWorld,
             parent,
             0,
             pset.get<bool>("foil_visible"),
             G4Colour::Blue(),
             pset.get<bool>("foil_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("GabsVacWinHolder",
             TubsParams(rVacWin, rIn, vacHolderThick/2),
             findMaterialOrThrow(pset.get<string>("wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, zend - lOut + vacHolderThick/2) - parent.centerInWorld,
             parent,
             0,
             pset.get<bool>("wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("GabsVacuumWindow",
             TubsParams(0., rVacWin, tVacWin/2),
             findMaterialOrThrow(pset.get<string>("vacuum_window_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, zend - exitFoilThick - foilDistance - tVacWin/2) - parent.centerInWorld,
             parent,
             0,
             pset.get<bool>("foil_visible"),
             G4Colour::Blue(),
             pset.get<bool>("foil_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    //----------------
    CLHEP::Hep3Vector  bulkGasRefPoint(0,0,0);
    vector<double> zPlanes = { zend - exitFoilThick - foilDistance,
                               zend - lOut + vacHolderThick,
                               zend - lOut + vacHolderThick,
                               zend - exitFoilThick};

    vector<double> r1{0., 0., 0., 0.};

    vector<double> r2{ rVacWin, rVacWin, rIn, rIn };

    nestPolycone("GabsBulkGas",
                 mu2e::PolyconsParams(zPlanes, r1, r2),
                 findMaterialOrThrow(pset.get<string>("inside_material")),
                 0, // no rotation
                 bulkGasRefPoint - parent.centerInWorld,
                 parent,
                 0,
                 pset.get<bool>("gas_visible"),
                 G4Colour::Yellow(),
                 pset.get<bool>("gas_solid"),
                 forceAuxEdgeVisible_,
                 placePV_,
                 doSurfaceCheck_
                 );


    
  }

  //================================================================
} // end namespace mucap
