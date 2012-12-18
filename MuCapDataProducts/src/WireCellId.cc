#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mu2e {

  WireCellId::WireCellId(const WirePlaneId& plane, unsigned int cell)
    : plane_(plane)
    , cell_(cell)
  {}


  std::ostream& operator<<(std::ostream& os, const WireCellId& id) {
    return os<<"WireCellId("<<id.plane().number()
             <<","<<id.cell()
             <<" )";
  }

  namespace { int MAGIC_FACTOR = 1000; } // less than 1000 wires in any plane

  unsigned WireCellId::encodeToInteger() const {
    return MAGIC_FACTOR * plane_.number() + cell_;
  }

  WireCellId decodeFromInteger(unsigned encodedVolumeId) {
    unsigned cell = encodedVolumeId % MAGIC_FACTOR;
    unsigned plane = encodedVolumeId / MAGIC_FACTOR;
    return WireCellId(WirePlaneId(plane), cell);
  }
}
