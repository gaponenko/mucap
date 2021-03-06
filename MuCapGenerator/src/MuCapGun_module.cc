// A generator for the cases when the desired distribution of particle
// properties can be factorized into the following independent
// distributions:
//
//      position,
//      momentum direction,
//      momentum magnitude
//      time
//
// Andrei Gaponenko, 2012

#include <string>
#include <memory>
#include <cmath>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

#include "SeedService/inc/SeedService.hh"

#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "MuCapInterfaces/inc/IPositionGenerator.hh"
#include "MuCapInterfaces/inc/IAngleGenerator.hh"
#include "MuCapInterfaces/inc/ISpectrumGenerator.hh"

#include "MuCapGenerator/inc/generatorHelpers.hh"

#include "TH1.h"

namespace mucap {
  using fhicl::ParameterSet;
  using mu2e::GenParticle;
  using mu2e::GenParticleCollection;
  using mu2e::PDGCode;

  class MuCapGun : public art::EDProducer {
    art::RandomNumberGenerator::base_engine_t& eng_;

    PDGCode::type pdgId_;
    double mass_;

    std::unique_ptr<IPositionGenerator> pos_;
    std::unique_ptr<ISpectrumGenerator> energyVar_;
    bool energyVarIsEk_;
    std::unique_ptr<IAngleGenerator> omega_;
    std::unique_ptr<ISpectrumGenerator> time_;

    TH1* hek_;
    TH1* hmomentum_;

  public:
    explicit MuCapGun(fhicl::ParameterSet const& pset);
    virtual void produce(art::Event& event);
  };

  MuCapGun::MuCapGun(fhicl::ParameterSet const& pset)
    : eng_(createEngine(art::ServiceHandle<mu2e::SeedService>()->getSeed()))
    , pdgId_(PDGCode::type(pset.get<int>("pdgId")))
    , mass_(mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , pos_(makePositionGenerator(pset.get<ParameterSet>("position"), eng_))
    , energyVarIsEk_(true)
    , omega_(makeAngleGenerator(pset.get<ParameterSet>("angles"), eng_))
    , time_(makeSpectrumGenerator(pset.get<ParameterSet>("time"), eng_))
    , hek_()
    , hmomentum_()
  {
    produces<GenParticleCollection>();

    const fhicl::ParameterSet energyPars = pset.get<ParameterSet>("energySpec");
    const std::string evname = energyPars.get<std::string>("variable");

    energyVar_ = makeSpectrumGenerator(energyPars, eng_);
    if(evname == "kineticEnergy") {
      energyVarIsEk_ = true;
    }
    else if(evname == "momentum") {
      energyVarIsEk_ = false;
    }
    else {
      throw cet::exception("CONFIG")<<__func__<<": unknown energySpec variable: "<<evname<<"\n";
    }

    // the parametrization is specifically for energy, momentum var would probably be due to a mistake
    if(energyPars.get<std::string>("spectrum") == "MECO") {
      if(!energyVarIsEk_) {
	throw cet::exception("CONFIG")<<__func__<<": unknown energySpec variable is not kineticEnergy for the MECO spectrum\n";
      }
    }

    ParameterSet hset(pset.get<ParameterSet>("histograms"));
    art::ServiceHandle<art::TFileService> tfs;

    {
      ParameterSet hpars(hset.get<ParameterSet>("kineticEnergy", ParameterSet()));
      if(!hpars.is_empty()) {
        hek_ = tfs->make<TH1D>("kineticEnergy", "kinetic energy",
                               hpars.get<unsigned>("nbins"),
                               hpars.get<double>("xmin"),
                               hpars.get<double>("xmax")
                               );
      }
    }

    {
      ParameterSet hpars(hset.get<ParameterSet>("momentum", ParameterSet()));
      if(!hpars.is_empty()) {
        hmomentum_ = tfs->make<TH1D>("momentum", "momentum",
                                     hpars.get<unsigned>("nbins"),
                                     hpars.get<double>("xmin"),
                                     hpars.get<double>("xmax")
                                     );
      }
    }

  }

  void MuCapGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    double ek=0., etot=0., pmag=0.;
    if(energyVarIsEk_) {
      ek = energyVar_->generate();
      etot = mass_ + ek;
      pmag = sqrt(std::pow(etot,2) - std::pow(mass_, 2));
    }
    else { // var is momentum
      pmag = energyVar_->generate();
      etot = std::sqrt(std::pow(pmag,2) + std::pow(mass_, 2));
      ek = etot - mass_;
    }
    if(hek_) { hek_->Fill(ek); }
    if(hmomentum_) { hmomentum_->Fill(pmag); }

    const CLHEP::HepLorentzVector momentum(pmag * omega_->generate(), etot);

    //std::cout<<"MuCapGun: generated p = "<<momentum<<std::endl;

    output->push_back(GenParticle(pdgId_,
                                  mu2e::GenId::particleGun,
                                  pos_->generate(),
                                  momentum,
                                  time_->generate()
                                  ));

    event.put(std::move(output));
  }

} // namespace mucap

DEFINE_ART_MODULE(mucap::MuCapGun);
