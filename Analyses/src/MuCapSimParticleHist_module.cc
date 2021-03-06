// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "TH1.h"
#include "TH2.h"

namespace mucap {

  using mu2e::SimParticle;
  using mu2e::SimParticleCollection;

  //================================================================
  class MuCapSimParticleHist : public art::EDAnalyzer {
    std::string inModuleLabel_;
    std::string inInstanceName_;

    TH1 *primary_ek_;
    TH1 *primary_mom_;
    TH1 *primary_z_;

    // Different ranges
    TH1 *interaction_z0_;
    TH1 *interaction_z1_;

    TH2 *primary_xy_;

  public:
    explicit MuCapSimParticleHist(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  MuCapSimParticleHist::MuCapSimParticleHist(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inInstanceName_(pset.get<std::string>("inputInstanceName", ""))
    , primary_ek_()
    , primary_mom_()
    , primary_z_()
    , interaction_z0_()
    , interaction_z1_()
    , primary_xy_()
  {}

  //================================================================
  void MuCapSimParticleHist::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

    primary_ek_ =  tfs->make<TH1D>("primary_ek", "Primary kinetic energy", 500, 0., 50.);
    primary_ek_->GetXaxis()->SetTitle("Ek [MeV]");

    primary_mom_ =  tfs->make<TH1D>("primary_p", "Primary momentum", 600, 0., 300.);
    primary_mom_->GetXaxis()->SetTitle("p [MeV/c]");

    primary_z_ =  tfs->make<TH1D>("primary_z", "Primary start z", 1000, -0.5, +0.5);
    primary_z_->GetXaxis()->SetTitle("z [mm]");

    interaction_z0_ =  tfs->make<TH1D>("interaction_z0", "Interaction z (zoomed)", 1000, -0.5, +0.5);
    interaction_z0_->GetXaxis()->SetTitle("z [mm]");

    interaction_z1_ =  tfs->make<TH1D>("interaction_z1", "Interaction z", 4800, -600, +600);
    interaction_z1_->GetXaxis()->SetTitle("z [mm]");

    primary_xy_ =  tfs->make<TH2D>("primary_xy", "Primary vtx y vs x", 201, -100.5, +100.5, 201, -100.5, +100.5);
    primary_xy_->GetXaxis()->SetTitle("x [mm]");
    primary_xy_->GetYaxis()->SetTitle("y [mm]");
    primary_xy_->SetOption("colz");
  }

  //================================================================
  void MuCapSimParticleHist::analyze(const art::Event& event) {

    art::Handle<SimParticleCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
    const SimParticleCollection& coll(*ih);

    if(!coll.empty()) {
      cet::map_vector_key iprim(1);
      const SimParticle primary(coll.getOrThrow(iprim));
      primary_mom_->Fill(primary.startMomentum().vect().mag());
      primary_ek_->Fill(primary.startMomentum().e() - primary.startMomentum().m());
      primary_z_->Fill(primary.startPosition().z());
      primary_xy_->Fill(primary.startPosition().x(), primary.startPosition().y());

      for(SimParticleCollection::const_iterator i=coll.begin(); i!=coll.end(); ++i) {
        if(i->first != iprim) {
          const SimParticle& sp(i->second);
          interaction_z0_->Fill(sp.startPosition().z());
          interaction_z1_->Fill(sp.startPosition().z());
        }
      }
    }
  }

  //================================================================
} // namespace mucap

DEFINE_ART_MODULE(mucap::MuCapSimParticleHist);
