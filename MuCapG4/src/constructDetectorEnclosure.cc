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
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

#include "MuCapG4/inc/MuCapSD.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mucap {

  using namespace std;
  using namespace mu2e;
  using fhicl::ParameterSet;

  mu2e::VolumeInfo MuCapWorld::constructDetectorEnclosure(const VolumeInfo& parent, const fhicl::ParameterSet& pset) {

    CLHEP::Hep3Vector detectorCenter(0,0,0);
    const double rInside = pset.get<double>("hous_radius") - pset.get<double>("hous_mantle");

    // The volume are arranged differently than in G3: for easier
    // visualization with the GDML viewer I don't want to make Al
    // walls a "container" volume.  Instead the "big" container volume
    // is filled with He/N2, and solid walls (leaf volumes) are placed
    // into it.

    TubsParams mother_params(0., pset.get<double>("hous_radius"), pset.get<double>("hous_length")/2+pset.get<double>("hous_cover"));

    VolumeInfo mother = nestTubs("HOUS_MOTHER",
                                 mother_params,
                                 findMaterialOrThrow(pset.get<string>("hous_inside_material")),
                                 0, // no rotation
                                 detectorCenter,
                                 parent,
                                 0,
                                 pset.get<bool>("hous_visible"),
                                 G4Colour::White(),
                                 pset.get<bool>("hous_solid"),
                                 forceAuxEdgeVisible_,
                                 placePV_,
                                 doSurfaceCheck_
                                 );

    TubsParams hous_wall_params(rInside,
                                pset.get<double>("hous_radius"),
                                pset.get<double>("hous_length")/2 + pset.get<double>("hous_cover"));

    VolumeInfo wall = nestTubs("HOUS_WALL",
                               hous_wall_params,
                               findMaterialOrThrow(pset.get<string>("hous_wall_material")),
                               0, // no rotation
                               CLHEP::Hep3Vector(0,0,0), // center in parent
                               mother,
                               0,
                               pset.get<bool>("hous_wall_visible"),
                               G4Colour::Grey(),
                               pset.get<bool>("hous_wall_solid"),
                               forceAuxEdgeVisible_,
                               placePV_,
                               doSurfaceCheck_
                               );

    //----------------------------------------------------------------
    // Internal (He/N2 side) part of the end cover

    const double hlen1 = pset.get<double>("hent_length_1");
    const double hent1_r = pset.get<double>("hent_radius_1");
    TubsParams hen1_ring_params(hent1_r, rInside, hlen1/2);
    const double hen1_ring_zabs = pset.get<double>("hous_length")/2 + hlen1/2;

    nestTubs("housInRingUp",
             hen1_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, -hen1_ring_zabs),
             mother,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("housInRingDn",
             hen1_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, +hen1_ring_zabs),
             mother,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    //----------------------------------------------------------------
    // Foils delimiting the He/N2 volume

    const double foil_thick = pset.get<double>("foil_thickness");
    TubsParams foil_params(0., hent1_r, foil_thick/2);
    const double foil_zabs = pset.get<double>("hous_length")/2 + hlen1 - foil_thick/2;

    nestTubs("housFoilUp",
             foil_params,
             findMaterialOrThrow(pset.get<string>("hous_foil_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, -foil_zabs),
             mother,
             0,
             pset.get<bool>("hous_foil_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_foil_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("housFoilDn",
             foil_params,
             findMaterialOrThrow(pset.get<string>("hous_foil_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, +foil_zabs),
             mother,
             0,
             pset.get<bool>("hous_foil_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_foil_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    //----------------------------------------------------------------
    // The part outside of the He/N2 volume.
    // First fill it with air, than place Al pieces inside air.

    const double hlen2 = pset.get<double>("hous_cover") - hlen1;
    TubsParams hen2_air_params(0, rInside, hlen2/2);
    const double hen2_air_zabs = pset.get<double>("hous_length")/2 + hlen1 + hlen2/2;

    VolumeInfo airUp = nestTubs("housAirUp",
                                hen2_air_params,
                                findMaterialOrThrow(pset.get<string>("hous_outside_material")),
                                0, // no rotation
                                CLHEP::Hep3Vector(0,0, -hen2_air_zabs),
                                mother,
                                0,
                                pset.get<bool>("hous_visible"),
                                G4Colour::Yellow(),
                                pset.get<bool>("hous_solid"),
                                forceAuxEdgeVisible_,
                                placePV_,
                                doSurfaceCheck_
                                );

    VolumeInfo airDn = nestTubs("housAirDn",
                                hen2_air_params,
                                findMaterialOrThrow(pset.get<string>("hous_outside_material")),
                                0, // no rotation
                                CLHEP::Hep3Vector(0,0, +hen2_air_zabs),
                                mother,
                                0,
                                pset.get<bool>("hous_visible"),
                                G4Colour::Yellow(),
                                pset.get<bool>("hous_solid"),
                                forceAuxEdgeVisible_,
                                placePV_,
                                doSurfaceCheck_
                                );

    //----------------
    // The air-side ring

    const double hent2_r = pset.get<double>("hent_radius_2");
    TubsParams hen2_ring_params(hent2_r, rInside, hlen2/2);

    nestTubs("housOutRingUp",
             hen2_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, 0),
             airUp,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("housOutRingDn",
             hen2_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, 0),
             airDn,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    //----------------
    // The retaining ring

    const double flunsh_len = pset.get<double>("hent_flunsh");
    TubsParams flunsh_ring_params(hent1_r, hent2_r, flunsh_len/2);
    const double flunsh_zabs = hlen2/2 - flunsh_len/2;

    nestTubs("housFlunshUp",
             flunsh_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, +flunsh_zabs),
             airUp,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    nestTubs("housFlunshDn",
             flunsh_ring_params,
             findMaterialOrThrow(pset.get<string>("hous_wall_material")),
             0, // no rotation
             CLHEP::Hep3Vector(0,0, -flunsh_zabs),
             airDn,
             0,
             pset.get<bool>("hous_wall_visible"),
             G4Colour::Grey(),
             pset.get<bool>("hous_wall_solid"),
             forceAuxEdgeVisible_,
             placePV_,
             doSurfaceCheck_
             );

    //----------------
    return mother;
  }

  //================================================================
} // end namespace mucap
