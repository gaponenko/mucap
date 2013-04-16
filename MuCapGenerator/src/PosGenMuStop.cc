// Andrei Gaponenko, 2013

#include "MuCapGenerator/inc/PosGenMuStop.hh"

#include <fstream>
#include <sstream>
#include <cmath>

#include "cetlib/exception.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MuCapGeom/inc/Geometry.hh"

#include "MuCapGeom/inc/Geometry.hh"

namespace mucap {

  PosGenMuStop::point PosGenMuStop::uv2xy(const point& uv) {
    static const double pi = 4*atan(1.);
    static const double rot = 45*pi/180.; // turns X axis into U axis
    static const double cosrot = cos(rot);
    static const double sinrot = sin(rot);

    return point(cosrot * uv.x - sinrot * uv.y,
                 sinrot * uv.x + cosrot * uv.y);
  }


  PosGenMuStop::PosGenMuStop(const std::string& inputFileName,
                             CLHEP::HepRandomEngine& rng)
    : flat_(rng)
    , zmin_()
    , zmax_()
  {
    art::ServiceHandle<Geometry> geom;
    zmin_ = geom->targetCenter().z() - 0.5*geom->targetThickness();
    zmax_ = geom->targetCenter().z() + 0.5*geom->targetThickness();

    const std::string muStopsFileName =  mu2e::ConfigFileLookupPolicy()(inputFileName);

    std::ifstream in(muStopsFileName.c_str());
    if(!in) {
      throw cet::exception("CONFIG")<<"PosGenMuStop: Can not open input file "<<muStopsFileName<<"\n";
    }

    std::string line;
    while(getline(in, line)) {
      //std::cout<<"Got line "<<line<<std::endl;
      std::istringstream is(line);
      double u=0, v=0;
      if(!(is>>u>>v)) {
        throw cet::exception("INPUTS")<<"Error: bad line in input file "<<muStopsFileName<<":\n"<<line<<std::endl;
      }
      xy_.push_back(uv2xy(point(u*CLHEP::centimeter, v*CLHEP::centimeter)));
    }

    if(xy_.empty()) {
      throw cet::exception("CONFIG")<<"PosGenMuStop: did not get any muon stops from "<<muStopsFileName<<"\n";
    }
    else {
      std::cout<<"PosGenMuStop: read "<<xy_.size()<<" muon stops from "<<muStopsFileName<<std::endl;
    }
  }

  //================================================================
  CLHEP::Hep3Vector PosGenMuStop::generate() {
    const double z = zmin_ + (zmax_ - zmin_)*flat_.fire();
    const int entry = flat_.fireInt(xy_.size());
    return CLHEP::Hep3Vector(xy_[entry].x, xy_[entry].y, z);
  }

}
