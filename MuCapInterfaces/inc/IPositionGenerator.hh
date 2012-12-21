// Andrei Gaponenko, 2012

#ifndef MuCapInterfaces_inc_IPositionGenerator_hh
#define MuCapInterfaces_inc_IPositionGenerator_hh

#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {
  class IPositionGenerator {
  public:

    // non-const because it modifies the random number generator state
    virtual CLHEP::Hep3Vector generate() = 0;

    ~IPositionGenerator() {}
  };
}

#endif/*MuCapInterfaces_inc_IPositionGenerator_hh*/
