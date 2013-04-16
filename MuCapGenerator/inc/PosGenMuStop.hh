// Draws transverse stopped positions from a sample of reconstructed
// muon stops in data.  The z position is uniform inside the stopping
// target.
//
// Andrei Gaponenko, 2013

#ifndef MuCapGenerator_inc_PosGenMuStop_hh
#define MuCapGenerator_inc_PosGenMuStop_hh

#include "MuCapInterfaces/inc/IPositionGenerator.hh"

#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {


  class PosGenMuStop: public IPositionGenerator {
  public:
    virtual CLHEP::Hep3Vector generate();

    PosGenMuStop(const std::string& muStopsFileName,
                 CLHEP::HepRandomEngine& rng);

  private:
    struct point { double x; double y; point(double a, double b) : x(a), y(b) {} };
    static point uv2xy(const point& uv);
    std::vector<point> xy_;
    CLHEP::RandFlat flat_;
    double zmin_;
    double zmax_;
  };
}

#endif/*MuCapGenerator_inc_PosGenMuStop_hh*/
