// Andrei Gaponenko, 2012

#ifndef MuCapGeom_inc_MuCapGeom_hh
#define MuCapGeom_inc_MuCapGeom_hh

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"

namespace mucap {
  class Geometry {
  public:
    Geometry(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~Geometry();

    double targetThickness() const;

  private:
  };
}

#endif/*MuCapGeom_inc_MuCapGeom_hh*/
