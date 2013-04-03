// Generate a spectrum using its tabulated values
//
// Andrei Gaponenko, 2013

#ifndef MuCapGenerator_inc_SpectrumGenTabulated_hh
#define MuCapGenerator_inc_SpectrumGenTabulated_hh

#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"
#include "CLHEP/Random/RandGeneral.h"


namespace mucap {

  class TabulatedFunction;

  class SpectrumGenTabulated: public ISpectrumGenerator {
  public:

    virtual double generate() override;

    // The function must have at least 2 points defined
    SpectrumGenTabulated(const TabulatedFunction& func, CLHEP::HepRandomEngine& rng);

  private:
    // The archaic CLHEP interface wants a double* to a table of values
    // as a constructor argument.  The otherwise unnecessary tmp_ here
    // is to keep the table while RandGeneral is being constructed.
    std::vector<double> tmp_;
    std::vector<double> getYValues(const TabulatedFunction& func);

    CLHEP::RandGeneral spectrum_;
    double xmin_;
    double xmax_;


  };
}

#endif/*MuCapGenerator_inc_SpectrumGenTabulated_hh*/
