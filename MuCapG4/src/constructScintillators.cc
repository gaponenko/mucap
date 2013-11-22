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

  void MuCapWorld::constructScintillators(const VolumeInfo& parent, const fhicl::ParameterSet& pset) {
    typedef std::vector<double> VD;
    const VD scz = pset.get<VD>("zcenter");
    const VD rIn = pset.get<VD>("rIn");
    if(rIn.size() != scz.size()) {
      throw cet::exception("GEOM")<<__func__<<" different lengths of input zcenter and rIn\n";
    }
    const VD rOut = pset.get<VD>("rOut");
    if(rOut.size() != scz.size()) {
      throw cet::exception("GEOM")<<__func__<<" different lengths of input zcenter and rOut\n";
    }
    const VD thick = pset.get<VD>("thick");
    if(thick.size() != scz.size()) {
      throw cet::exception("GEOM")<<__func__<<" different lengths of input zcenter and thick\n";
    }

    const VD twrap = pset.get<VD>("wrap_thicknesses");
    const vector<string> mwrap = pset.get<vector<string> >("wrap_materials");
    if(twrap.size() != mwrap.size()) {
      throw cet::exception("GEOM")<<__func__<<" different lengths of input wrap_thicknesses and wrap_materials\n";
    }

    for(unsigned isc = 0; isc < scz.size(); ++isc) {
      TubsParams scint_params(rIn[isc], rOut[isc], thick[isc]/2);
      std::ostringstream os;
      os<<"scintillator_"<<isc;

      CLHEP::Hep3Vector scintCenter = CLHEP::Hep3Vector(0, 0, scz[isc]) - parent.centerInWorld;

      nestTubs(os.str(),
               scint_params,
               findMaterialOrThrow(pset.get<string>("scint_material")),
               0, // no rotation
               scintCenter,
               parent,
               0,
               pset.get<bool>("scint_visible"),
               G4Colour::Magenta(),
               pset.get<bool>("scint_solid"),
               forceAuxEdgeVisible_,
               placePV_,
               doSurfaceCheck_
               );

      double alreadyWrappedThickness = 0.;
      for(unsigned iwrap = 0; iwrap < twrap.size(); ++iwrap) {
        std::ostringstream ow;
        ow<<"scint_"<<isc<<"_wrap_"<<iwrap;
        TubsParams wrap_params(rIn[isc], rOut[isc], twrap[iwrap]/2);
        const CLHEP::Hep3Vector
          wrapOffset(0, 0,
                     scint_params.zHalfLength() + alreadyWrappedThickness + wrap_params.zHalfLength());

        nestTubs(ow.str()+"_up",
                 wrap_params,
                 findMaterialOrThrow(mwrap[iwrap]),
                 0, // no rotation
                 scintCenter - wrapOffset,
                 parent,
                 0,
                 pset.get<bool>("wrap_visible"),
                 G4Colour::Black(),
                 pset.get<bool>("wrap_solid"),
                 forceAuxEdgeVisible_,
                 placePV_,
                 doSurfaceCheck_
                 );

        nestTubs(ow.str()+"_dn",
                 wrap_params,
                 findMaterialOrThrow(mwrap[iwrap]),
                 0, // no rotation
                 scintCenter + wrapOffset,
                 parent,
                 0,
                 pset.get<bool>("wrap_visible"),
                 G4Colour::Black(),
                 pset.get<bool>("wrap_solid"),
                 forceAuxEdgeVisible_,
                 placePV_,
                 doSurfaceCheck_
                 );
        
        alreadyWrappedThickness += twrap[iwrap];
      }
    }
  }

  //================================================================
} // end namespace mucap
