#include "MuCapDataProducts/inc/MuCapRecoHit.hh"

namespace mucap {
  std::ostream& operator<<(std::ostream& os, const MuCapRecoHit& hit) {
    return os<<"MuCapRecoHit( "<<hit.rid()<<", t = "<<hit.time()<<", w = "<<hit.width()<<" )";
  }
}
