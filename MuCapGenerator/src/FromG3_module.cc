// Re-use proton kinematics from a G3 run for G3-G4 comparisons.
//
// Andrei Gaponenko, 2013

#include <string>
#include <memory>
#include <cmath>
#include <fstream>

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

#include "TH1.h"

namespace mucap {
  using fhicl::ParameterSet;
  using mu2e::GenParticle;
  using mu2e::GenParticleCollection;
  using mu2e::PDGCode;

  namespace {
    struct G3Info {
      CLHEP::Hep3Vector startPos;
      CLHEP::Hep3Vector startMom;
      CLHEP::Hep3Vector endPos;
    };

    bool getline(std::istream& inp, G3Info *res) {
      std::string buf;
      if(std::getline(inp,buf)) {
      std::istringstream is(buf);
        double x_start(0), y_start(0), z_start(0), p_x(0), p_y(0), p_z(0), x_stop(0), y_stop(0), z_stop(0);
        if(is>>x_start>>y_start>>z_start>>p_x>>p_y>>p_z>>x_stop>>y_stop>>z_stop) {
          res->startPos = CLHEP::Hep3Vector(x_start, y_start, z_start);
          res->startMom = CLHEP::Hep3Vector(p_x, p_y, p_z);
          res->endPos = CLHEP::Hep3Vector(x_stop, y_stop, z_stop);
          // Convert cm to mm
          res->startPos *= 10.;
          res->endPos *= 10.;
          // Convert GeV/c to MeV/c
          res->startMom *= 1000.;
          return true;
        }
      }
      return false;
    }
  }

  class FromG3 : public art::EDProducer {
    std::string inputFileName_;
    std::ifstream inp_;

    PDGCode::type pdgId_;
    double mass_;

    TH1* hek_;
    TH1* hmomentum_;

  public:
    explicit FromG3(fhicl::ParameterSet const& pset);
    virtual void produce(art::Event& event);
  };

  FromG3::FromG3(fhicl::ParameterSet const& pset)
    : inputFileName_(pset.get<std::string>("inputFileName"))
    , inp_(inputFileName_.c_str())
    , pdgId_(PDGCode::type(pset.get<int>("pdgId")))
    , mass_(mu2e::GlobalConstantsHandle<mu2e::ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , hek_()
    , hmomentum_()
  {
    produces<GenParticleCollection>();
    produces<GenParticleCollection>("G3ProtonStops");

    if(!inp_.is_open()) {
      throw cet::exception("CONFIG")<<"Can not open input file "<<inputFileName_<<"\n";
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

  void FromG3::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    std::unique_ptr<GenParticleCollection> g3stops(new GenParticleCollection);

    G3Info buf;
    if(getline(inp_, &buf)) {

      const double pmag = buf.startMom.mag();
      if(hmomentum_) { hmomentum_->Fill(pmag); }
      const double etot = sqrt(std::pow(mass_, 2) + std::pow(pmag, 2));
      const double ek = etot - mass_;
      if(hek_) { hek_->Fill(ek); }

      const CLHEP::HepLorentzVector momentum(buf.startMom, etot);

      output->push_back(GenParticle(pdgId_,
                                    mu2e::GenId::particleGun,
                                    buf.startPos,
                                    momentum,
                                    0
                                    ));

      // We just want to store the stop position from G3, but it's easiest to
      // put it in GenParticle
      g3stops->push_back(GenParticle(pdgId_,
                                     mu2e::GenId::particleGun,
                                     buf.endPos,
                                     CLHEP::HepLorentzVector(),
                                     0
                                     ));
    }

    event.put(std::move(output));
    event.put(std::move(g3stops), "G3ProtonStops");
  }

} // namespace mucap

DEFINE_ART_MODULE(mucap::FromG3);
