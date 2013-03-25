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

    std::auto_ptr<IPositionGenerator> pos_;
    std::auto_ptr<ISpectrumGenerator> ek_;
    std::auto_ptr<IAngleGenerator> omega_;
    std::auto_ptr<ISpectrumGenerator> time_;

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
    , ek_(makeSpectrumGenerator(pset.get<ParameterSet>("kineticEnergy"), eng_))
    , omega_(makeAngleGenerator(pset.get<ParameterSet>("angles"), eng_))
    , time_(makeSpectrumGenerator(pset.get<ParameterSet>("time"), eng_))
    , hek_()
    , hmomentum_()
  {
    produces<GenParticleCollection>();

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
    std::auto_ptr<GenParticleCollection> output(new GenParticleCollection);

    const double ek = ek_->generate();
    if(hek_) { hek_->Fill(ek); }

    const double etot = mass_ + ek;
    const double pmag = sqrt(std::pow(etot,2) - std::pow(mass_, 2));
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
