#include "MuCapGeom/inc/Geometry.hh"

#include "art/Framework/Services/Registry/ServiceMacros.h"

#include <string>
#include <vector>

namespace mucap {

  using namespace std;
  using namespace fhicl;

  namespace {
    const string TARGET_MODULE_TYPE = "targetModule";
  }

  Geometry::Geometry(fhicl::ParameterSet const& pset, art::ActivityRegistry&iRegistry)
    : pset_(pset)
    , targetCenter_()
    , targetThickness_(pset.get<ParameterSet>(TARGET_MODULE_TYPE)
                       .get<ParameterSet>("target")
                       .get<double>("thickness"))
  {
    //----------------------------------------------------------------
    vector<ParameterSet> modulePars(pset_.get<vector<ParameterSet> >("chamberModules"));
    unsigned itgt(-1);
    for(unsigned i = 0; i < modulePars.size(); ++i) {
      if(modulePars[i].get<string>("type") == TARGET_MODULE_TYPE) {
        if(itgt == -1u) {
          itgt = i;
        }
        else {
          throw cet::exception("GEOM")<<"Multiple modules of type "<<TARGET_MODULE_TYPE<<" are not allowed\n";
        }
      }
    }

    if(itgt == -1u) {
      throw cet::exception("GEOM")<<"A module of type "<<TARGET_MODULE_TYPE<<" is required\n";
    }

    const vector<double> zfoil(modulePars[itgt].get<vector<double> >("zfoil"));
    if(zfoil.size()%2 != 1) {
      throw cet::exception("GEOM")<<__func__<<" module of type "<<TARGET_MODULE_TYPE
                                  <<" must contain an odd number of foils.  Got: "<<zfoil.size()
                                  <<"\n";
    }
    const unsigned iTargetFoil = zfoil.size()/2;
    targetCenter_ = CLHEP::Hep3Vector(0, 0, zfoil[iTargetFoil]);

  } // end Geometry() ctr

} // end namespace mucap

DEFINE_ART_SERVICE(mucap::Geometry);
