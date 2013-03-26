#include "MuCapDataProducts/inc/WireReadoutId.hh"

namespace mucap {

  WireReadoutId::WireReadoutId(const WirePlaneId& plane, unsigned int channel)
    : plane_(plane)
    , channel_(channel)
  {}

  std::ostream& operator<<(std::ostream& os, const WireReadoutId& id) {
    return os<<"WireReadoutId( "<<id.plane().number()<<", "<<id.channel()<<" )";
  }

}
