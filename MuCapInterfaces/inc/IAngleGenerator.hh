// Andrei Gaponenko, 2012

#ifndef MuCapInterfaces_inc_IAngleGenerator_hh
#define MuCapInterfaces_inc_IAngleGenerator_hh

#include "CLHEP/Vector/ThreeVector.h"

namespace mucap {
  class IAngleGenerator {
  public:

    // Non-const because it modifies the random number generator state
    // The output is a unit vector.
    virtual CLHEP::Hep3Vector generate() = 0;

    virtual ~IAngleGenerator() {}
  };
}

#endif/*MuCapInterfaces_inc_IAngleGenerator_hh*/
