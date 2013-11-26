#include "MuCapGenerator/inc/generatorHelpers.hh"

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <limits>

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h" // for twopi

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "MuCapGeom/inc/Geometry.hh"

#include "MuCapGenerator/inc/PosGenCylinder.hh"
#include "MuCapGenerator/inc/PosGenMuStop.hh"
#include "MuCapGenerator/inc/AngleGenUniform.hh"
#include "MuCapGenerator/inc/SpectrumGenFlat.hh"
#include "MuCapGenerator/inc/SpectrumGenMECO.hh"
#include "MuCapGenerator/inc/SpectrumGenExp.hh"
#include "MuCapGenerator/inc/SpectrumGenGauss.hh"
#include "MuCapGenerator/inc/SpectrumGenTabulated.hh"

#include "MuCapUtilities/inc/mecoSpectrum.hh"
#include "MuCapUtilities/inc/TabulatedFunction.hh"


namespace mucap {

  //================================================================
  std::unique_ptr<IPositionGenerator>
  makePositionGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {

    const std::string shape(pset.get<std::string>("shape"));

    //----------------------------------------------------------------
    if(shape == "cylinder") {
      return std::unique_ptr<IPositionGenerator>
        (new PosGenCylinder(pset.get<CLHEP::Hep3Vector>("position"),
                            pset.get<double>("radius"),
                            pset.get<double>("halfdz"),
                            eng));
    }
    //----------------------------------------------------------------
    else if(shape == "targetCylinder") {

      art::ServiceHandle<Geometry> geom;

      return std::unique_ptr<IPositionGenerator>
        (new PosGenCylinder(geom->targetCenter(),
                            pset.get<double>("radius"),
                            0.5*geom->targetThickness(),
                            eng));
    }
    //----------------------------------------------------------------
    else if(shape == "muonStop") {
      return std::unique_ptr<IPositionGenerator>
        (new PosGenMuStop(pset.get<std::string>("fileName"), eng));
    }

    //----------------------------------------------------------------
    throw cet::exception("GEOM")<<__func__<<": unknown shape setting \""<<shape<<"\"\n";

  } // makePositionGenerator()

  //================================================================
  std::unique_ptr<ISpectrumGenerator>
  makeSpectrumGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {

    const std::string spectrum(pset.get<std::string>("spectrum"));

    //----------------------------------------------------------------
    if(spectrum == "flat") {
      return std::unique_ptr<ISpectrumGenerator>
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
      pars.T0inv = 1./pset.get<double>("T0");

      return std::unique_ptr<ISpectrumGenerator>
        (new SpectrumGenMECO(pars, eng, pset.get<unsigned>("maxIter", 1000)));
    }
    //----------------------------------------------------------------
    else if(spectrum == "exp") {
      return std::unique_ptr<ISpectrumGenerator>
        (new SpectrumGenExp(pset.get<double>("scale"),
                            pset.get<double>("min", 0),
                            pset.get<double>("max", std::numeric_limits<double>::max()),
                            eng));
    }
    //----------------------------------------------------------------
    else if(spectrum == "gauss") {
      return std::unique_ptr<ISpectrumGenerator>
        (new SpectrumGenGauss(pset.get<double>("mean"),
                              pset.get<double>("sigma"),
                              eng));
    }
    //----------------------------------------------------------------
    else if(spectrum == "tabulated") {
      mu2e::ConfigFileLookupPolicy mu2eFileFinder;
      TabulatedFunction func(mu2eFileFinder(pset.get<std::string>("tableFileName")));
      return std::unique_ptr<ISpectrumGenerator>(new SpectrumGenTabulated(func, eng));
    }
    //----------------------------------------------------------------
    throw cet::exception("GEOM")<<__func__<<": unknown spectrum setting \""<<spectrum<<"\"\n";

  } // makeSpectrumGenerator()

  //================================================================
  std::unique_ptr<IAngleGenerator>
  makeAngleGenerator(const fhicl::ParameterSet& pset, art::RandomNumberGenerator::base_engine_t& eng) {
    return std::unique_ptr<IAngleGenerator>
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
