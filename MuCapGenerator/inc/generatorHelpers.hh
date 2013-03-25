// Andrei Gaponenko, 2012

#ifndef MuCapGenerator_inc_generatorHelpers_hh
#define MuCapGenerator_inc_generatorHelpers_hh

#include <memory>

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace fhicl { class ParameterSet; }

namespace mucap {

  //================================================================
  class IPositionGenerator;

  std::unique_ptr<IPositionGenerator>
  makePositionGenerator(const fhicl::ParameterSet& pset,
                        art::RandomNumberGenerator::base_engine_t& eng);

  //================================================================
  class ISpectrumGenerator;

  std::unique_ptr<ISpectrumGenerator>
  makeSpectrumGenerator(const fhicl::ParameterSet& pset,
                        art::RandomNumberGenerator::base_engine_t& eng);

  //================================================================
  class IAngleGenerator;

  std::unique_ptr<IAngleGenerator>
  makeAngleGenerator(const fhicl::ParameterSet& pset,
                        art::RandomNumberGenerator::base_engine_t& eng);

  //================================================================

}

#endif/*MuCapGenerator_inc_generatorHelpers_hh*/
