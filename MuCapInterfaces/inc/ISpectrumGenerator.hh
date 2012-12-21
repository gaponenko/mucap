// Distribution of any one-dimensional variable.
//
// Andrei Gaponenko, 2012

#ifndef MuCapInterfaces_inc_ISpectrumGenerator_hh
#define MuCapInterfaces_inc_ISpectrumGenerator_hh

namespace mucap {
  class ISpectrumGenerator {
  public:

    // non-const because it modifies the random number generator state
    virtual double generate() = 0;

    ~ISpectrumGenerator() {}
  };
}

#endif/*MuCapInterfaces_inc_ISpectrumGenerator_hh*/
