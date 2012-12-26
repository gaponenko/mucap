#include "MuCapGenerator/inc/generatorHelpers.hh"

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h" // for twopi

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

#include "MuCapGeom/inc/Geometry.hh"

#include "MuCapGenerator/inc/PosGenCylinder.hh"
#include "MuCapGenerator/inc/AngleGenUniform.hh"
#include "MuCapGenerator/inc/SpectrumGenFlat.hh"
#include "MuCapGenerator/inc/SpectrumGenMECO.hh"

#include "MuCapUtilities/inc/mecoSpectrum.hh"


namespace mucap {

  //================================================================
  std::auto_ptr<IPositionGenerator>
  makePositionGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {

    const std::string shape(pset.get<std::string>("shape"));

    //----------------------------------------------------------------
    if(shape == "cylinder") {
      return std::auto_ptr<IPositionGenerator>
        (new PosGenCylinder(pset.get<CLHEP::Hep3Vector>("position"),
                            pset.get<double>("radius"),
                            pset.get<double>("halfdz"),
                            eng));
    }
    //----------------------------------------------------------------
    else if(shape == "targetCylinder") {

      art::ServiceHandle<Geometry> geom;

      return std::auto_ptr<IPositionGenerator>
        (new PosGenCylinder(geom->targetCenter(),
                            pset.get<double>("radius"),
                            0.5*geom->targetThickness(),
                            eng));
    }

    //----------------------------------------------------------------
    throw cet::exception("GEOM")<<__func__<<": unknown shape setting \""<<shape<<"\"\n";

  } // makePositionGenerator()

  //================================================================
  std::auto_ptr<ISpectrumGenerator>
  makeSpectrumGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {

    const std::string spectrum(pset.get<std::string>("spectrum"));

    //----------------------------------------------------------------
    if(spectrum == "flat") {
      return std::auto_ptr<ISpectrumGenerator>
        (new SpectrumGenFlat(pset.get<double>("center"),
                             pset.get<double>("halfWidth"),
                             eng));
    }
    //----------------------------------------------------------------
    else if(spectrum == "MECO") {

      MECOPars pars;
      pars.A = 1;
      pars.Tth = pset.get<double>("Tth");
      pars.alpha = pset.get<double>("alpha");
      pars.T0 = pset.get<double>("T0");

      return std::auto_ptr<ISpectrumGenerator>
        (new SpectrumGenMECO(pars, eng, pset.get<unsigned>("maxIter", 1000)));
    }
    //----------------------------------------------------------------
    throw cet::exception("GEOM")<<__func__<<": unknown spectrum setting \""<<spectrum<<"\"\n";

  } // makeSpectrumGenerator()

  //================================================================
  std::auto_ptr<IAngleGenerator>
  makeAngleGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {
    return std::auto_ptr<IAngleGenerator>
      (new AngleGenUniform(mu2e::RandomUnitSphereParams(pset.get<double>("czmin"),
                                                        pset.get<double>("czmax"),
                                                        pset.get<double>("phimin", 0.),
                                                        pset.get<double>("phimax", CLHEP::twopi)
                                                        )
                           , eng)
       );

  } // makeAngleGenerator()

  //================================================================
}
