// Andrei Gaponenko, 2013

#include "cetlib/exception.h"

// Mu2e includes.
#include "MuCapDataProducts/inc/WirePlaneId.hh"
#include "MuCapDataProducts/inc/WireReadoutId.hh"
#include "MuCapDataProducts/inc/MuCapRecoHit.hh"
#include "MuCapDataProducts/inc/MuCapRecoHitCollection.hh"

#include "MuCapGeom/inc/Geometry.hh"

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

  const unsigned MAX_PLANES_D = 44;
  const unsigned MAX_WIRES_D = 80;
  const unsigned MAX_PLANES_P = 12;
  const unsigned MAX_WIRES_P = 160;
  const unsigned MAXHITS_DC = MAX_WIRES_D * MAX_PLANES_D * 2; // dc_mngw
  const unsigned MAXHITS_PC = MAX_WIRES_P * MAX_PLANES_P * 2; // pc_mngw

  // Why do DetectorGeo.h defs disagree with MOFIA's det_geom_mod.f90??? Use here constants from the latter source.
  const unsigned AG_MAXM_SCINTS = 8;
  const unsigned AG_MAXT_SCINTS = 8;
  const unsigned AG_MAXO_SCINTS = 64;
  const unsigned AG_MAX_SCINTS = AG_MAXM_SCINTS + AG_MAXT_SCINTS + AG_MAXO_SCINTS;
  const unsigned MAXHITS_SC = AG_MAX_SCINTS * 8; // sc_mngw

  //================================================================
  // EVID branch
  struct Clark_EVID {
    int nrun;
    int nevt;
    Clark_EVID() : nrun(0), nevt(0) {}
  };

  //================================================================
  // Event branch
  struct Clark_Event {

    Int_t treeversion;
    Int_t timestamp;
    Int_t type;
    Float_t m12width;
    Float_t cptime[3];
    Float_t rftime[3];
    Int_t nwin;
    Int_t ntr;
    Int_t pienuitr;

    Clark_Event()
      : treeversion(0)
      , timestamp(0)
      , type(0)
      , m12width(0)
      , nwin(0)
      , ntr(0)
      , pienuitr(0)
    {
      cptime[0]=cptime[1]=cptime[2]=0.;
      rftime[0]=rftime[1]=rftime[2]=0.;
    }
  };

  //================================================================
  struct Clark_WireHits {
  };

  //================================================================
  class ClarkTreeDumper : public art::EDAnalyzer {
    fhicl::ParameterSet pset_;
    std::string digiModuleLabel_;
    art::ServiceHandle<Geometry> geom_;

    // Tree branch variables

    Clark_EVID evid_;
    Clark_Event event_;

    //_________________________ MuCapture: All hits - DC  __________________________//
    Int_t dc_nhits_;
    Float_t dc_time_[MAXHITS_DC];//[dc_nhits]
    Float_t dc_width_[MAXHITS_DC];//[dc_nhits]
    Int_t dc_plane_[MAXHITS_DC];//[dc_nhits]
    Int_t dc_cell_[MAXHITS_DC];//[dc_nhits]

    //_________________________ MuCapture: All hits - PC  __________________________//
    Int_t pc_nhits_;
    Float_t pc_time_[MAXHITS_PC];//[pc_nhits]
    Float_t pc_width_[MAXHITS_PC];//[pc_nhits]
    Int_t pc_plane_[MAXHITS_PC];//[pc_nhits]
    Int_t pc_cell_[MAXHITS_PC];//[pc_nhits]

    //_________________________ MuCapture: All hits - SC  __________________________//
    Int_t sc_nhits_;
    Float_t sc_time_[MAXHITS_SC];//[sc_nhits]
    Float_t sc_width_[MAXHITS_SC];//[sc_nhits]
    Int_t sc_iscint_[MAXHITS_SC];//[sc_nhits]
    Int_t sc_wire_[MAXHITS_SC];//[sc_nhits]

    //----------------
    TTree *nt_;

  public:
    explicit ClarkTreeDumper(const fhicl::ParameterSet& pset);
    virtual void beginJob() override;
    virtual void analyze(const art::Event& event) override;
  };

  //================================================================
  ClarkTreeDumper::ClarkTreeDumper(const fhicl::ParameterSet& pset)
    : pset_(pset)
    , digiModuleLabel_(pset.get<std::string>("digiModuleLabel"))
    , dc_nhits_(0)
    , pc_nhits_(0)
    , sc_nhits_(0)
    , nt_(nullptr)
  {
  }

  //================================================================
  void ClarkTreeDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    nt_ = tfs->make<TTree>( "T", "Hits for Clark");
    nt_->Branch("EVID", &evid_, "nrun/I:nevt/I");
    nt_->Branch("Event", &event_, "treeversion/I:timestamp/I:type/I:m12width/F:cptime[3]/F:rftime[3]/F:nwin/I:ntr/I:pienuitr/I");

    TBranch *bdchits = nt_->Branch("DC_hits", 0, "DC_nhits/I:DC_time[DC_nhits]/F:DC_width[DC_nhits]/F:DC_plane[DC_nhits]/I:DC_cell[DC_nhits]/I");
    bdchits->FindLeaf("DC_nhits")->SetAddress(&dc_nhits_);
    bdchits->FindLeaf("DC_time")->SetAddress(&dc_time_);
    bdchits->FindLeaf("DC_width")->SetAddress(&dc_width_);
    bdchits->FindLeaf("DC_plane")->SetAddress(&dc_plane_);
    bdchits->FindLeaf("DC_cell")->SetAddress(&dc_cell_);

    TBranch *bpchits = nt_->Branch("PC_hits", 0, "PC_nhits/I:PC_time[PC_nhits]/F:PC_width[PC_nhits]/F:PC_plane[PC_nhits]/I:PC_cell[PC_nhits]/I");
    bpchits->FindLeaf("PC_nhits")->SetAddress(&pc_nhits_);
    bpchits->FindLeaf("PC_time")->SetAddress(&pc_time_);
    bpchits->FindLeaf("PC_width")->SetAddress(&pc_width_);
    bpchits->FindLeaf("PC_plane")->SetAddress(&pc_plane_);
    bpchits->FindLeaf("PC_cell")->SetAddress(&pc_cell_);

    TBranch *bschits = nt_->Branch("SC_hits", 0, "SC_nhits/I:SC_time[SC_nhits]/F:SC_width[SC_nhits]/F:SC_iscint[SC_nhits]/I:SC_wire[SC_nhits]/I");
    bschits->FindLeaf("SC_nhits")->SetAddress(&sc_nhits_);
    bschits->FindLeaf("SC_time")->SetAddress(&sc_time_);
    bschits->FindLeaf("SC_width")->SetAddress(&sc_width_);
    bschits->FindLeaf("SC_iscint")->SetAddress(&sc_iscint_);
    bschits->FindLeaf("SC_wire")->SetAddress(&sc_wire_);
  }

  //================================================================
  void ClarkTreeDumper::analyze(const art::Event& event) {
    // The mofia tree uses (run,event) EVID.  There is no room to store the 3-component art EventID.
    // We store subRun, and drop the run number.
    evid_.nrun = event.subRun();
    evid_.nevt = event.event();

    //----------------------------------------------------------------
    art::Handle<MuCapRecoHitCollection> hdchits;
    event.getByLabel(digiModuleLabel_, "dchits", hdchits);
    const MuCapRecoHitCollection& dchits(*hdchits);
    if(dchits.size() > MAXHITS_DC) {
      throw cet::exception("OVERFLOW") << "ClarkTreeDumper: too many DC hits\n";
    }

    dc_nhits_ = dchits.size();
    for(unsigned i=0; i<dchits.size(); ++i) {
      dc_time_[i] = dchits[i].time();
      dc_width_[i] = dchits[i].width();
      const WireReadoutId rid = dchits[i].rid();
      dc_plane_[i] = geom_->wpByGlobalNumber(rid.plane().number()).localPlane;
      dc_cell_[i] = rid.channel();
    }

    //----------------------------------------------------------------
    art::Handle<MuCapRecoHitCollection> hpchits;
    event.getByLabel(digiModuleLabel_, "pchits", hpchits);
    const MuCapRecoHitCollection& pchits(*hpchits);
    if(pchits.size() > MAXHITS_PC) {
      throw cet::exception("OVERFLOW") << "ClarkTreeDumper: too many PC hits\n";
    }

    pc_nhits_ = pchits.size();
    for(unsigned i=0; i<pchits.size(); ++i) {
      pc_time_[i] = pchits[i].time();
      pc_width_[i] = pchits[i].width();
      const WireReadoutId rid = pchits[i].rid();
      pc_plane_[i] = geom_->wpByGlobalNumber(rid.plane().number()).localPlane;
      pc_cell_[i] = rid.channel();
    }

    //----------------------------------------------------------------
    nt_->Fill();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mucap::ClarkTreeDumper);
