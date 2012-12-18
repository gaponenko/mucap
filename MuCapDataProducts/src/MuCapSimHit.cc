#include "MuCapDataProducts/inc/MuCapSimHit.hh"

#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  MuCapSimHit::MuCapSimHit(const StepPointMC& sp)
    : sp_(&sp)
    , cid_()
  {}

  std::ostream& operator<<(std::ostream& os, const MuCapSimHit& sh) {
    return os<<"MuCapSimHit("<<sh.cid()<<", edep"<<sh.hit().totalEDep()<<" )";
  }
}
