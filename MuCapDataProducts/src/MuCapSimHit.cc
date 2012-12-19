#include "MuCapDataProducts/inc/MuCapSimHit.hh"

#include "MCDataProducts/inc/StepPointMC.hh"

namespace mucap {

  MuCapSimHit::MuCapSimHit(const mu2e::StepPointMC& sp)
    : sp_(&sp)
    , cid_(WireCellId::decodeFromInteger(sp.volumeId()))
  {}

  std::ostream& operator<<(std::ostream& os, const MuCapSimHit& sh) {
    return os<<"MuCapSimHit( "<<sh.cid()<<", edep = "<<sh.hit().totalEDep()<<" )";
  }
}
