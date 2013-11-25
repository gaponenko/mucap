// Andrei Gaponenko, 2013

#include "MuCapG4/inc/MuCapWorld.hh"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <utility>

#include "cetlib/exception.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Paraboloid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
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
#include "G4PVPlacement.hh"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include "MuCapG4/inc/MuCapSD.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mucap {

  using namespace std;
  using namespace mu2e;
  using fhicl::ParameterSet;
  using mu2e::VolumeInfo;
  using CLHEP::Hep3Vector;

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
    constructTEC(vac, bpipe_params.innerRadius(), pset.get<ParameterSet>("tec"));
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


    const double bulge = pset.get<double>("vacuum_window_bulge");
    if(bulge < 0.) {
      throw cet::exception("GEOM")<<__func__<<" gabs.vacuum_window_bulge < 0\n";
    }
    else if(bulge > 0.) {
      //----------------------------------------------------------------
      // Bulged window and gas

      CLHEP::Hep3Vector gasCenter(0,0, zend - exitFoilThick - foilDistance - bulge/2);
      VolumeInfo bulgeGas("GabsBulgeGas", gasCenter - parent.centerInWorld, parent.centerInWorld);
      bulgeGas.solid = new G4Paraboloid(bulgeGas.name, bulge/2, 0., rVacWin);

      finishNesting(bulgeGas,
                    findMaterialOrThrow(pset.get<string>("inside_material")),
                    0, // no rotation
                    bulgeGas.centerInParent,
                    parent.logical,
                    0,
                    pset.get<bool>("gas_visible"),
                    G4Colour::Yellow(),
                    pset.get<bool>("gas_solid"),
                    forceAuxEdgeVisible_,
                    placePV_,
                    doSurfaceCheck_
                    );

      // Info for the bulged window
      CLHEP::Hep3Vector winCenter(0,0, zend - exitFoilThick - foilDistance - (bulge + tVacWin)/2);
      VolumeInfo vacuumWindow("GabsVacuumWindow", winCenter - parent.centerInWorld, parent.centerInWorld);

      // A taller (than the gas bulge) paraboloid following the same curve as the bulge.
      G4Paraboloid* vwenvelope = new G4Paraboloid("vwenvelope", (bulge + tVacWin)/2, 0., rVacWin*sqrt(1 + tVacWin/bulge));

      // This is the bulged foil, but it unwanted margins at radius > rVacWin
      G4SubtractionSolid *vwsub = new G4SubtractionSolid("vwsub", vwenvelope, bulgeGas.solid, 0, Hep3Vector(0, 0, tVacWin/2));

      // A cylinder to cut off the radial margins.  The precise length does not matter, just make it long enough.
      G4Tubs* vwcut = new G4Tubs("vwcut", 0., rVacWin, tVacWin + bulge, 0., CLHEP::twopi);

      vacuumWindow.solid = new G4IntersectionSolid("GabsVacuumWindow", vwsub, vwcut, 0, Hep3Vector(0,0,0));
      finishNesting(vacuumWindow,
                    findMaterialOrThrow(pset.get<string>("vacuum_window_material")),
                    0, // no rotation
                    vacuumWindow.centerInParent,
                    parent.logical,
                    0,
                    pset.get<bool>("foil_visible"),
                    G4Colour::Blue(),
                    pset.get<bool>("foil_solid"),
                    forceAuxEdgeVisible_,
                    placePV_,
                    doSurfaceCheck_
                    );
    }
    else {
      //----------------------------------------------------------------
      // Flat window

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
    }

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
  void MuCapWorld::constructTEC(const VolumeInfo& parent, double rOut, const fhicl::ParameterSet& pset) {
    if(pset.get<bool>("installed")) {
      const double halfLength = pset.get<double>("length")/2;

      G4Material *gasMaterial = findMaterialOrThrow(pset.get<string>("gas_material"));
      const VolumeInfo gas = nestTubs("TECGas",
                                      TubsParams(0., rOut, halfLength),
                                      gasMaterial,
                                      0, // no rotation
                                      Hep3Vector(0,0, pset.get<double>("center_z")) - parent.centerInWorld,
                                      parent,
                                      0,
                                      pset.get<bool>("gas_visible"),
                                      G4Colour::Yellow(),
                                      pset.get<bool>("gas_solid"),
                                      forceAuxEdgeVisible_,
                                      placePV_,
                                      doSurfaceCheck_
                                      );

      //----------------
      const double foilHalfThick(pset.get<double>("foil_thickness")/2);
      TubsParams foilParams(0., rOut, foilHalfThick);

      nestTubs("TECFoilUp",
               foilParams,
               findMaterialOrThrow(pset.get<string>("foil_material")),
               0, // no rotation
               Hep3Vector(0,0, -(halfLength-foilHalfThick)),
               gas,
               0,
               pset.get<bool>("foil_visible"),
               G4Colour::Yellow(),
               pset.get<bool>("foil_solid"),
               forceAuxEdgeVisible_,
               placePV_,
               doSurfaceCheck_
               );

      nestTubs("TECFoilDn",
               foilParams,
               findMaterialOrThrow(pset.get<string>("foil_material")),
               0, // no rotation
               Hep3Vector(0,0, +(halfLength-foilHalfThick)),
               gas,
               0,
               pset.get<bool>("foil_visible"),
               G4Colour::Yellow(),
               pset.get<bool>("foil_solid"),
               forceAuxEdgeVisible_,
               placePV_,
               doSurfaceCheck_
               );

      //----------------------------------------------------------------
      // Place field wire planes
      const vector<double> fwz = pset.get<vector<double> >("field_plane_zoffset");
      const vector<double> fwr = pset.get<vector<double> >("field_plane_rotation");
      if(fwz.size() != fwr.size()) {
        throw cet::exception("GEOM")<<__func__
                                    <<": field_plane_zoffset and field_plane_rotation must be of the same size\n";
      }

      G4LogicalVolume* fieldPlane = makeTECFieldPlane(rOut, gasMaterial, pset);
      for(unsigned i=0; i<fwz.size(); ++i) {
        std::ostringstream os;
        os<<"TECFieldPlane"<<i;
        new G4PVPlacement(HepGeom::RotateZ3D(fwr[i]*CLHEP::degree)*HepGeom::TranslateZ3D(fwz[i]),
                          fieldPlane,
                          os.str(),
                          gas.logical,
                          false,
                          i,
                          doSurfaceCheck_);
      }

      //----------------------------------------------------------------

    }
  }

  //================================================================
  G4LogicalVolume *MuCapWorld::makeTECFieldPlane(double radius, G4Material *gasMaterial, const fhicl::ParameterSet& pset) {
    const double wireRadius(pset.get<double>("field_wire_diameter")/2);
    G4LogicalVolume *plane = new G4LogicalVolume(new G4Tubs("TECFieldPlane", 0, radius, wireRadius, 0, CLHEP::twopi),
                                                 gasMaterial,
                                                 "TECFieldPlane"
                                                 );


    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    CLHEP::HepRotation *wireRotationInv = reg.add(new CLHEP::HepRotation());
    wireRotationInv->rotateY(90*CLHEP::degree);
    
    G4Material *bronze = findMaterialOrThrow(pset.get<string>("field_wire_material"));
    const double pitch = pset.get<double>("field_wire_pitch");
    double offset =0;
    for(int i=0; (offset = (0.5+i)*pitch) + wireRadius < radius; ++i) {
      std::ostringstream os;
      os<<"TECFieldWire"<<i;
      const double halfLength = sqrt(std::pow(radius,2) - std::pow(offset + wireRadius, 2));

      nestTubs(os.str()+"a",
               TubsParams(0, wireRadius, halfLength),
               bronze,
               wireRotationInv,
               Hep3Vector(0,offset,0),
               plane,
               0,
               pset.get<bool>("field_wire_visible"),
               G4Colour::Red(),
               pset.get<bool>("field_wire_solid"),
               forceAuxEdgeVisible_,
               placePV_,
               doSurfaceCheck_
               );

      nestTubs(os.str()+"b",
               TubsParams(0, wireRadius, halfLength),
               bronze,
               wireRotationInv,
               Hep3Vector(0,-offset,0),
               plane,
               0,
               pset.get<bool>("field_wire_visible"),
               G4Colour::Red(),
               pset.get<bool>("field_wire_solid"),
               forceAuxEdgeVisible_,
               placePV_,
               doSurfaceCheck_
               );
    }

    return plane;
  }

  //================================================================
} // end namespace mucap
