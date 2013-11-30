// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"

#include "MuCapDataProducts/inc/WirePlaneId.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"
#include "MuCapDataProducts/inc/MuCapSimHit.hh"
#include "MuCapDataProducts/inc/WireReadoutId.hh"
#include "MuCapDataProducts/inc/MuCapRecoHit.hh"
#include "MuCapDataProducts/inc/MuCapRecoHitCollection.hh"

#include "MuCapDigitization/inc/MuCapCharge.hh"

#include "MuCapGeom/inc/Geometry.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mucap {

  //================================================================

  using mu2e::StepPointMCCollection;

  class MuCapDigitization : public art::EDProducer {
    fhicl::ParameterSet pset_;
    std::string hitsModuleLabel_;
    std::string hitsInstanceName_;
    std::string particlesModuleLabel_;
    std::string particlesInstanceName_;
    double minTimeWidthDC_;
    double minTimeWidthPC_;
    double daqGateTimeMin_;
    double daqGateTimeMax_;

    art::ServiceHandle<Geometry> geom_;

    void digitize(MuCapRecoHitCollection *dchits, MuCapRecoHitCollection *pchits, const StepPointMCCollection& inputs);
    WireReadoutId getReadoutId(const WireCellId& cid);
    
    TH1 *hnumSimHitsPerChannelDC_;
    TH1 *hsimHitTimeSepDC_;
    TH1 *hnumSimHitsPerChannelPC_;
    TH1 *hsimHitTimeSepPC_;

  public:
    explicit MuCapDigitization(const fhicl::ParameterSet& pset);
    virtual void beginJob() override;
    virtual void produce(art::Event& event) override;
  };

  //================================================================
  MuCapDigitization::MuCapDigitization(const fhicl::ParameterSet& pset)
    : pset_(pset)
    , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
    , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))
    , particlesModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
    , particlesInstanceName_(pset.get<std::string>("particlesInstanceName", ""))

    , minTimeWidthDC_(pset.get<double>("minHitTimeSepDC"))
    , minTimeWidthPC_(pset.get<double>("minHitTimeSepPC"))
    , daqGateTimeMin_(pset.get<double>("daqGateTimeMin"))
    , daqGateTimeMax_(pset.get<double>("daqGateTimeMax"))

    , hnumSimHitsPerChannelDC_(nullptr)
    , hsimHitTimeSepDC_(nullptr)
    , hnumSimHitsPerChannelPC_(nullptr)
    , hsimHitTimeSepPC_(nullptr)
  {
    produces<MuCapRecoHitCollection>("dchits");
    produces<MuCapRecoHitCollection>("pchits");
  }

  //================================================================
  void MuCapDigitization::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    hnumSimHitsPerChannelDC_ =  tfs->make<TH1D>("numSimHitsPerChannelDC", "numSimHitsPerChannel, DC", 50, -0.5, 49.5);
    hsimHitTimeSepDC_ =  tfs->make<TH1D>("simHitTimeSepDC", "simHitTimeSep, DC", 400, 0., 2000.);
    hnumSimHitsPerChannelPC_ =  tfs->make<TH1D>("numSimHitsPerChannelPC", "numSimHitsPerChannel, PC", 50, -0.5, 49.5);
    hsimHitTimeSepPC_ =  tfs->make<TH1D>("simHitTimeSepPC", "simHitTimeSep, PC", 400, 0., 2000.);
  }

  //================================================================
  void MuCapDigitization::produce(art::Event& event) {
    art::Handle<StepPointMCCollection> inputs;
    event.getByLabel(hitsModuleLabel_, hitsInstanceName_, inputs);

    std::unique_ptr<MuCapRecoHitCollection> dchits(new MuCapRecoHitCollection());
    std::unique_ptr<MuCapRecoHitCollection> pchits(new MuCapRecoHitCollection());

    digitize(&*dchits, &*pchits, *inputs);

    event.put(std::move(dchits), "dchits");
    event.put(std::move(pchits), "pchits");
  }

  //================================================================
  void MuCapDigitization::digitize(MuCapRecoHitCollection *dchits,
                                   MuCapRecoHitCollection *pchits,
                                   const StepPointMCCollection& inputs)
  {
    // Sort sim hits by readout channel and time
    MuCapChargeCollection cc;
    for(unsigned i=0; i<inputs.size(); ++i) {
      const MuCapSimHit sh(inputs[i]);
      const WireReadoutId rid(getReadoutId(sh.cid()));
      const double registeredTime = 0/*FIXME: drift*/ + sh.hit().time();
      cc[rid].push(MuCapTimedChargeDeposit(registeredTime, sh.hit().ionizingEdep()));
    }

    // Merge hits on the same channel that are close in time, and form output hits
    for(auto& iter : cc) {
      const WireReadoutId& rid = iter.first;
      const unsigned globalPlane = rid.plane().number();
      const WPType wpt = geom_->wpByGlobalNumber(globalPlane).wpt;
      TH1* hnum = (wpt == WPType::DC) ?  hnumSimHitsPerChannelDC_ : hnumSimHitsPerChannelPC_;
      TH1* hsep = (wpt == WPType::DC) ?  hsimHitTimeSepDC_ : hsimHitTimeSepPC_;
      const double tcut = (wpt == WPType::DC) ? minTimeWidthDC_: minTimeWidthPC_;
      MuCapRecoHitCollection *out = (wpt == WPType::DC) ? dchits : pchits;

      MuCapChargeHistory& hq = iter.second;
      while(!hq.empty()) {
        // Start a new hit
        const double tstart = hq.top().time;
        double t = tstart;
        double esum = hq.top().energy;
        int numSimHits = 1;

        hq.pop();
        while(!hq.empty() && (hq.top().time - t <= tcut)) {
          hsep->Fill(hq.top().time - t);
          t = hq.top().time;
          esum += hq.top().energy;
          ++numSimHits;
          hq.pop();
        }
        hnum->Fill(numSimHits);

        if((daqGateTimeMin_ <= tstart) && (tstart < daqGateTimeMax_)) {
          out->emplace_back(rid, tstart, (t-tstart)+tcut);
        }
      }
    }
  }

  //================================================================
  WireReadoutId MuCapDigitization::getReadoutId(const WireCellId& cid) {
    // FIXME: need to impelement ganged wires
    return WireReadoutId(cid.plane(), cid.cell());
  }

  //================================================================
} // namespace mucap

DEFINE_ART_MODULE(mucap::MuCapDigitization);
