// Andrei Gaponenko, 2012, based on Mu2eStudyWorld by Krzysztof Genser.

#include "MuCapG4/inc/MuCapWorld.hh"

// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/constructStudyEnv_v001.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
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
  MuCapWorld::MuCapWorld(const fhicl::ParameterSet& geom)
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

    const fhicl::ParameterSet world(geom_.get<fhicl::ParameterSet>("world"));

    // the canonical World volume
    VolumeInfo worldVInfo(nestBox("World",
                                  world.get<vector<double> >("halfLength"),
                                  findMaterialOrThrow(world.get<string>("materialName")),
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
    vector<fhicl::ParameterSet> modulePars(geom_.get<vector<fhicl::ParameterSet> >("chamberModules"));
    for(unsigned i = 0, wpoff=0; i < modulePars.size(); ++i) {
      wpoff += constructChamberModule(i, wpoff, modulePars[i], worldVInfo);
    }

    //----------------------------------------------------------------

    return worldVInfo.physical;
  }

  //================================================================
  unsigned MuCapWorld::constructChamberModule(unsigned moduleNumber,
                                              unsigned planeNumberOffset,
                                              const fhicl::ParameterSet& pars,
                                              const VolumeInfo& parent
                                              )
  {
    std::cout<<"MuCapWorld: constructing module "<<moduleNumber
             <<", planeNumberOffset = "<<planeNumberOffset<<std::endl;

    const string moduleType(pars.get<string>("type"));
    const fhicl::ParameterSet detail(geom_.get<fhicl::ParameterSet>(moduleType));
    const fhicl::ParameterSet cathode(geom_.get<fhicl::ParameterSet>("cathode"));

    const fhicl::ParameterSet
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

    ostringstream osname("chamberModule_");
    osname<<std::setw(2)<<std::setfill('0')<<moduleNumber;
    VolumeInfo modInfo(nestTubs(osname.str(),
                                params,
                                findMaterialOrThrow(detail.get<string>("material")),
                                0, // no rotation
                                moduleCenterInParent,
                                parent,
                                0,
                                detail.get<bool>("visible"),
                                G4Colour::Cyan(),
                                detail.get<bool>("solid"),
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

    return zwire.size();
  }

  //================================================================
  void MuCapWorld::constructFoil(unsigned imodule,
                                 unsigned ifoil,
                                 const fhicl::ParameterSet& foilPars,
                                 const CLHEP::Hep3Vector& centerInParent,
                                 const VolumeInfo& parent)
  {
  }

  //================================================================

} // end namespace mu2e
