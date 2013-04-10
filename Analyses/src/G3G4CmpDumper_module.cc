// Andrei Gaponenko, 2013

#include "cetlib/exception.h"

// Mu2e includes.
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticle.hh"

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/FindOne.h"
#include "art/Utilities/InputTag.h"

// ROOT
#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

// externals
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

// std
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mucap {

  using mu2e::SimParticleCollection;
  using mu2e::SimParticle;
  using mu2e::GenParticleCollection;
  using mu2e::GenParticle;
  
  //================================================================
  struct G3G4CmpInfo {
    double x_start{0.}, y_start{0.}, z_start{0.};
    double p_x{0.}, p_y{0.}, p_z{0.};
    double g4_xstop{0.}, g4_ystop{0.}, g4_zstop{0.};
    double g3_xstop{0.}, g3_ystop{0.}, g3_zstop{0.};
  };

  //================================================================
  class G3G4CmpDumper : public art::EDAnalyzer {
    std::string g4ModuleLabel_;
    std::string g4InstanceName_;
    std::string g3StopsModuleLabel_;
    std::string g3StopsInstanceName_;

    // Tree branch variables
    G3G4CmpInfo buf_;
    TTree *nt_;

  public:
    explicit G3G4CmpDumper(const fhicl::ParameterSet& pset);
    virtual void beginJob() override;
    virtual void analyze(const art::Event& event) override;
  };

  //================================================================
  G3G4CmpDumper::G3G4CmpDumper(const fhicl::ParameterSet& pset)
    : g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel"))
    , g4InstanceName_(pset.get<std::string>("g4InstanceName"))
    , g3StopsModuleLabel_(pset.get<std::string>("g3StopsModuleLabel"))
    , g3StopsInstanceName_(pset.get<std::string>("g3StopsInstanceName"))
    , nt_(nullptr)
  {}

  //================================================================
  void G3G4CmpDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    nt_ = tfs->make<TTree>( "cmp", "G3-G4 cmp");
    nt_->Branch("b", &buf_, "x_start/D:y_start/D:z_start/D:p_x/D:p_y/D:p_z/D:g4_xstop/D:g4_ystop/D:g4_zstop/D:g3_xstop/D:g3_ystop/D:g3_zstop/D");
  }

  //================================================================
  void G3G4CmpDumper::analyze(const art::Event& event) {
    art::Handle<SimParticleCollection> hg4particles;
    event.getByLabel(g4ModuleLabel_, g4InstanceName_, hg4particles);
    const SimParticleCollection& g4particles(*hg4particles);

    art::Handle<GenParticleCollection> g3particles;
    event.getByLabel(g3StopsModuleLabel_, g3StopsInstanceName_, g3particles);

    if(!g4particles.empty()) {
      if(!g3particles->size() == 1) {
        throw cet::exception("INPUTS")<<"Unexpected size of g3stops collection\n";
      }

      cet::map_vector_key iprim(1);
      const SimParticle primary(g4particles.getOrThrow(iprim));
      buf_.x_start = primary.startPosition().x();
      buf_.y_start = primary.startPosition().y();
      buf_.z_start = primary.startPosition().z();

      buf_.p_x = primary.startMomentum().x();
      buf_.p_y = primary.startMomentum().y();
      buf_.p_z = primary.startMomentum().z();

      buf_.g4_xstop = primary.endPosition().x();
      buf_.g4_ystop = primary.endPosition().y();
      buf_.g4_zstop = primary.endPosition().z();

      const GenParticle& g3part(g3particles->front());
      buf_.g3_xstop = g3part.position().x();
      buf_.g3_ystop = g3part.position().y();
      buf_.g3_zstop = g3part.position().z();

      nt_->Fill();
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mucap::G3G4CmpDumper);
