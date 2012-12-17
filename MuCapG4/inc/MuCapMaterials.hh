// Andrei Gaponenko, 2012

#ifndef MuCapG4_inc_MuCapMaterials_hh
#define MuCapG4_inc_MuCapMaterials_hh

#include "fhiclcpp/ParameterSet.h"

class G4Element;
class G4String;

namespace mu2e {

  class MuCapMaterials {
  public:

    explicit MuCapMaterials(const fhicl::ParameterSet& materialDefs) :
      pset_(materialDefs)
    {}

    // Construct all of the materials.
    // This is the interface of the class used by WorldMaker.
    void construct();

  private:
    fhicl::ParameterSet pset_;

    G4Element* getElementOrThrow(const G4String& name);
  };

} // end namespace mu2e
#endif /* MuCapG4_inc_MuCapMaterials_hh */
