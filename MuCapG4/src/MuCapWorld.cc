// Andrei Gaponenko, 2012, based on Mu2eStudyWorld by Krzysztof Genser.

#include "MuCapG4/inc/MuCapWorld.hh"

// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <utility>

// Framework includes
#include "cetlib/exception.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/constructStudyEnv_v001.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

// G4 includes
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
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

#include "Mu2eG4/inc/Mu2eGlobalField.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)


namespace mu2e {

  using namespace std;
  using fhicl::ParameterSet;

  const std::string MuCapWorld::targetModuleType_("targetModule");

  //================================================================
  MuCapWorld::MuCapWorld(const ParameterSet& geom)
    : geom_(geom)
    , forceAuxEdgeVisible_(geom.get<bool>("forceAuxEdgeVisible"))
    , doSurfaceCheck_(geom.get<bool>("doSurfaceCheck"))
    , placePV_(true)
  {}

  //================================================================
  MuCapWorld::~MuCapWorld(){

    // Do not destruct the solids, logical volumes or physical volumes.
    // G4 looks after that itself.

  }

  //================================================================
  // This is the callback called by G4
  G4VPhysicalVolume * MuCapWorld::construct(){

    const ParameterSet world(geom_.get<ParameterSet>("world"));

    // the canonical World volume
    VolumeInfo worldVInfo(nestBox("World",
                                  world.get<vector<double> >("halfLength"),
                                  findMaterialOrThrow(world.get<string>("material")),
                                  0,
                                  G4ThreeVector(),
                                  0, // no parent
                                  0, // we assign this volume copy number 0
                                  world.get<bool>("visible"),
                                  G4Colour::Cyan(),
                                  world.get<bool>("solid"),
                                  forceAuxEdgeVisible_,
                                  placePV_,
                                  false)); // do not surface check this one

    //----------------------------------------------------------------
    vector<ParameterSet> modulePars(geom_.get<vector<ParameterSet> >("chamberModules"));
    for(unsigned i = 0, wpoff=0; i < modulePars.size(); ++i) {
      wpoff += constructChamberModule(i, wpoff, modulePars[i], worldVInfo);
    }

    //----------------------------------------------------------------

    return worldVInfo.physical;
  }

  //================================================================
  //Helper for constructChamberModule()
  namespace {
    std::pair<const ParameterSet*, const ParameterSet*>
    getFoilPair(unsigned iplane, const ParameterSet& cathode, const ParameterSet& centralFoil, unsigned iTargetFoil) {
      return make_pair( (iplane==iTargetFoil)? &centralFoil : &cathode,
                        (iplane+1==iTargetFoil)? &centralFoil : &cathode
                        );
    }
  }

  //----------------
  unsigned MuCapWorld::constructChamberModule(unsigned moduleNumber,
                                              unsigned planeNumberOffset,
                                              const ParameterSet& pars,
                                              const VolumeInfo& parent
                                              )
  {
    std::cout<<"MuCapWorld: constructing module "<<moduleNumber
             <<", planeNumberOffset = "<<planeNumberOffset<<std::endl;

    const string moduleType(pars.get<string>("type"));
    const ParameterSet detail(geom_.get<ParameterSet>(moduleType));
    const ParameterSet cathode(geom_.get<ParameterSet>("cathode"));

    const ParameterSet
      centralFoil( (moduleType == targetModuleType_) ?
                   detail.get<ParameterSet>("target"):
                   cathode
                   );

    const vector<double> zfoil(pars.get<vector<double> >("zfoil"));
    const vector<double> zwire(pars.get<vector<double> >("zwire"));

    if(zwire.empty()) {
      throw cet::exception("GEOM")<<__func__<<" module with 0 wire planes requested\n";
    }

    if(zfoil.size() != 1+zwire.size()) {
      throw cet::exception("GEOM")<<__func__<<" inconsistent number of foils "
                                  <<zfoil.size()<<" and wire planes "<<zwire.size()
                                  <<"\n";
    }

    const double zcenter = 0.5*(zfoil.back() + zfoil.front());
    const double halfdz  = 0.5*(zfoil.back() - zfoil.front() + cathode.get<double>("thickness"));

    // module R == cathode R
    TubsParams params(0./* rIn */, cathode.get<double>("radius"), halfdz);

    const CLHEP::Hep3Vector moduleCenterInParent(0,0,zcenter);

    ostringstream osname;
    osname<<"chamberModule_"<<std::setw(2)<<std::setfill('0')<<moduleNumber;

    VolumeInfo modInfo(nestTubs(osname.str(),
                                params,
                                findMaterialOrThrow(detail.get<string>("material")),
                                0, // no rotation
                                moduleCenterInParent,
                                parent,
                                0,
                                detail.get<bool>("moduleBoxVisible"),
                                G4Colour::Cyan(),
                                detail.get<bool>("moduleBoxSolid"),
                                forceAuxEdgeVisible_,
                                placePV_,
                                doSurfaceCheck_
                                ));

    // Install the foils
    if((moduleType == targetModuleType_) && !(zfoil.size()%2)) {
      throw cet::exception("GEOM")<<__func__<<" module of type "<<targetModuleType_
                                  <<" must contain an odd number of foils.  Got: "<<zfoil.size()
                                  <<"\n";
    }
    const unsigned iTargetFoil = zfoil.size()/2;

    for(unsigned ifoil = 0; ifoil < zfoil.size(); ++ifoil) {
      const CLHEP::Hep3Vector foilCenterInParent(0,0, zfoil[ifoil] - zcenter);
      constructFoil(moduleNumber, ifoil,
                    ((ifoil==iTargetFoil)? centralFoil : cathode),
                    foilCenterInParent,
                    modInfo
                    );
    }

    // Create the drift cells and install wires
    const vector<double> planeRotationDegrees(pars.get<vector<double> >("rotation"));

    for(unsigned iplane = 0; iplane < zwire.size(); ++iplane) {

      std::pair<const ParameterSet*, const ParameterSet*>
        foils = getFoilPair(iplane, cathode, centralFoil, iTargetFoil);

      const double driftZmin = zfoil[iplane]   + 0.5*foils.first->get<double>("thickness") - zcenter;
      const double driftZmax = zfoil[iplane+1] - 0.5*foils.second->get<double>("thickness") - zcenter;

      constructDriftPlane(planeNumberOffset+iplane,
                          detail,
                          planeRotationDegrees[iplane] * CLHEP::degree,
                          zwire[iplane] - zcenter,
                          driftZmin,
                          driftZmax,
                          modInfo
                          );
    }



    return zwire.size();
  }

  //================================================================
  void MuCapWorld::constructFoil(unsigned imodule,
                                 unsigned ifoil,
                                 const ParameterSet& foilPars,
                                 const CLHEP::Hep3Vector& centerInParent,
                                 const VolumeInfo& parent)
  {
    ostringstream osname;
    osname<<"cathode_"<<std::setw(2)<<std::setfill('0')<<imodule<<"_"<<ifoil;

    TubsParams params(0./* rIn */, foilPars.get<double>("radius"), 0.5*foilPars.get<double>("thickness"));

    VolumeInfo modInfo(nestTubs(osname.str(),
                                params,
                                findMaterialOrThrow(foilPars.get<string>("material")),
                                0, // no rotation
                                centerInParent,
                                parent,
                                0,
                                foilPars.get<bool>("visible"),
                                G4Colour::Gray(),
                                foilPars.get<bool>("solid"),
                                forceAuxEdgeVisible_,
                                placePV_,
                                doSurfaceCheck_
                                ));
  }

  //================================================================
  void MuCapWorld::constructDriftPlane(unsigned globalPlaneNumber,
                                       const ParameterSet& detail,
                                       double planeRotationAngle,
                                       double zwire,
                                       double driftZmin,
                                       double driftZmax,
                                       const VolumeInfo& parent
                                       )
    {
      const unsigned ncells = detail.get<unsigned>("nwires");
      const double dx = detail.get<double>("wireSpacing");

      const double chamberRadius = geom_.get<ParameterSet>("cathode").get<double>("radius");
      const double safety = 0.5; // mm

      using namespace CLHEP;
      AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
      CLHEP::HepRotation *rot = reg.add(new HepRotation(HepRotationZ(-planeRotationAngle)));
      CLHEP::HepRotation *wireInCellRotation = reg.add(new HepRotation(HepRotationX(90*CLHEP::degree)));

      for(unsigned icell=0; icell<ncells; ++icell) {
        const double xmin = dx*(icell -0.5*ncells);
        const double xmax = xmin + dx;
        const double yHalfSize = -safety +
          sqrt(chamberRadius*chamberRadius - max(pow(xmin,2), pow(xmax,2)));

        if(yHalfSize <= 0) {
          throw cet::exception("GEOM")<<__func__<<" Got invalid yHalfSize = "<<yHalfSize<<"\n";
        }

        // std::cout<<"drift cell: x=["<<xmin<<", "<<xmax<<"], yHalfSize="<<yHalfSize<<std::endl;

        vector<double> cellHalfSize(3);
        cellHalfSize[0] = 0.5*dx;
        cellHalfSize[1] = yHalfSize;
        cellHalfSize[2] = 0.5*(driftZmax - driftZmin);

       const CLHEP::Hep3Vector
         cellCenterInParent(rot->inverse()*CLHEP::Hep3Vector(0.5*(xmax+xmin), 0, 0.5*(driftZmax + driftZmin)));

        ostringstream osname;
        osname<<"driftPlane_"<<std::setw(2)<<std::setfill('0')<<globalPlaneNumber<<"_"<<icell;
        VolumeInfo cellVI(nestBox(osname.str(),
                                  cellHalfSize,
                                  findMaterialOrThrow(detail.get<string>("material")),
                                  rot,
                                  cellCenterInParent,
                                  parent,
                                  1000*globalPlaneNumber + icell, // volume copy number
                                  detail.get<bool>("driftCellVisible"),
                                  G4Colour::Cyan(),
                                  detail.get<bool>("driftCellSolid"),
                                  forceAuxEdgeVisible_,
                                  placePV_,
                                  doSurfaceCheck_
                                  )
                          );

        // FIXME
        //cellVI.centerInWorld = cellCenterInParent + parent.centerInWorld;

        // Install the wire
        const ParameterSet wire(detail.get<ParameterSet>("wire"));
        TubsParams wirePars(0./* rIn */, 0.5*wire.get<double>("diameter"), yHalfSize);
        ostringstream wirename;
        wirename<<"wire_"<<std::setw(5)<<std::setfill('0')<<1000*globalPlaneNumber + icell;
        nestTubs(wirename.str(),
                 wirePars,
                 findMaterialOrThrow(wire.get<string>("material")),
                 wireInCellRotation,
                 Hep3Vector(0, 0, zwire - cellCenterInParent.z()),
                 cellVI,
                 1000*globalPlaneNumber + icell, // copy number
                 wire.get<bool>("visible"),
                 G4Colour::Yellow(),
                 wire.get<bool>("solid"),
                 forceAuxEdgeVisible_,
                 placePV_,
                 doSurfaceCheck_
                 );

      } //  for(cells)
    }

  //================================================================
} // end namespace mu2e
