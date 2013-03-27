#include "MuCapGeom/inc/Geometry.hh"


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

    , wpi_{
      { WPType::PC, 1},
      { WPType::PC, 2},
      { WPType::PC, 3},
      { WPType::PC, 4},

      { WPType::DC, 1},
      { WPType::DC, 2},
      { WPType::DC, 3},
      { WPType::DC, 4},
      { WPType::DC, 5},
      { WPType::DC, 6},
      { WPType::DC, 7},
      { WPType::DC, 8},
      { WPType::DC, 9},
      { WPType::DC, 10},
      { WPType::DC, 11},
      { WPType::DC, 12},
      { WPType::DC, 13},
      { WPType::DC, 14},
      { WPType::DC, 15},
      { WPType::DC, 16},
      { WPType::DC, 17},
      { WPType::DC, 18},
      { WPType::DC, 19},
      { WPType::DC, 20},
      { WPType::DC, 21},
      { WPType::DC, 22},

      { WPType::PC, 5},
      { WPType::PC, 6},
      { WPType::PC, 7},
      { WPType::PC, 8},

      { WPType::DC, 23},
      { WPType::DC, 24},
      { WPType::DC, 25},
      { WPType::DC, 26},
      { WPType::DC, 27},
      { WPType::DC, 28},
      { WPType::DC, 29},
      { WPType::DC, 30},
      { WPType::DC, 31},
      { WPType::DC, 32},
      { WPType::DC, 33},
      { WPType::DC, 34},
      { WPType::DC, 35},
      { WPType::DC, 36},
      { WPType::DC, 37},
      { WPType::DC, 38},
      { WPType::DC, 39},
      { WPType::DC, 40},
      { WPType::DC, 41},
      { WPType::DC, 42},
      { WPType::DC, 43},
      { WPType::DC, 44},

      { WPType::PC, 9},
      { WPType::PC, 10},
      { WPType::PC, 11},
      { WPType::PC, 12},
    }

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
