// Andrei Gaponenko, 2013

#include "cetlib/exception.h"

// Mu2e includes.
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MuCapDataProducts/inc/WirePlaneId.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"
#include "MuCapDataProducts/inc/MuCapSimHit.hh"

#include "GeometryService/inc/GeomHandle.hh"

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

  using mu2e::StepPointMCCollection;

  //================================================================
  struct ProtonHits {
    double pmag;
    double costh;
    double phi;

    int numHitPlanes;
    int hitPlaneMin;
    int hitPlaneMax;

    ProtonHits()
      : pmag(std::numeric_limits<double>::quiet_NaN())
      , costh(std::numeric_limits<double>::quiet_NaN())
      , phi(std::numeric_limits<double>::quiet_NaN())
      , numHitPlanes(-1)
      , hitPlaneMin(-1)
      , hitPlaneMax(-1)
    {}
  };


  //================================================================
  class ProtonHitDumper : public art::EDAnalyzer {
    fhicl::ParameterSet pset_;
    std::string hitsModuleLabel_;
    std::string hitsInstanceName_;
    std::string particlesModuleLabel_;
    std::string particlesInstanceName_;

    ProtonHits buf_;
    TTree *nt_;

  public:
    explicit ProtonHitDumper(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  ProtonHitDumper::ProtonHitDumper(const fhicl::ParameterSet& pset)
    : pset_(pset)
    , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
    , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))
    , particlesModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
    , particlesInstanceName_(pset.get<std::string>("particlesInstanceName", ""))
    , nt_(0)
  {}


  //================================================================
  void ProtonHitDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    static const char branchDesc[] = "pmag/D:costh/D:phi/D:numHitPlanes/I:hitPlaneMin/I:hitPlaneMax/I";
    nt_ = tfs->make<TTree>( "tr", "Proton hits ntuple");
    nt_->Branch("evt", &buf_, branchDesc);
  }

  //================================================================
  void ProtonHitDumper::analyze(const art::Event& event) {

    art::Handle<StepPointMCCollection> ih;
    event.getByLabel(hitsModuleLabel_, hitsInstanceName_, ih);
    const StepPointMCCollection& coll(*ih);

    if(!coll.empty()) {

      art::Handle<mu2e::SimParticleCollection> ih;
      event.getByLabel(particlesModuleLabel_, particlesInstanceName_, ih);

      cet::map_vector_key iprim(1);
      const mu2e::SimParticle primary(ih->getOrThrow(iprim));

      CLHEP::Hep3Vector mom = primary.startMomentum();
      buf_.pmag = mom.mag();
      buf_.costh = mom.cosTheta();
      buf_.phi = mom.phi();

      typedef std::set<WirePlaneId> PlaneSet;
      PlaneSet hitPlanes;

      for(unsigned i=0; i<coll.size(); ++i) {
        const MuCapSimHit sh(coll[i]);
        hitPlanes.insert(sh.cid().plane());
      }

      buf_.numHitPlanes = hitPlanes.size();
      buf_.hitPlaneMin = hitPlanes.begin()->number();
      buf_.hitPlaneMax = hitPlanes.rbegin()->number();

      nt_->Fill();
    }

  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mucap::ProtonHitDumper);
