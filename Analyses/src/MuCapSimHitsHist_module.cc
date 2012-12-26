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

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MuCapDataProducts/inc/WirePlaneId.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"
#include "MuCapDataProducts/inc/MuCapSimHit.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace mucap {

  using mu2e::StepPointMCCollection;

  //================================================================
  class MuCapSimHitsHist : public art::EDAnalyzer {
    fhicl::ParameterSet pset_;
    std::string hitsModuleLabel_;
    std::string hitsInstanceName_;
    std::string particlesModuleLabel_;
    std::string particlesInstanceName_;

    TH1 *cellEDepTotal_;
    TH1 *cellEDepIonizing_;
    TH2 *numCellPlane2d_;
    TH3 *numCellPlane3d_;

    double getPrimaryEk(const art::Event& event) const;

  public:
    explicit MuCapSimHitsHist(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  MuCapSimHitsHist::MuCapSimHitsHist(const fhicl::ParameterSet& pset)
    : pset_(pset)
    , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
    , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))
    , particlesModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
    , particlesInstanceName_(pset.get<std::string>("particlesInstanceName", ""))
    , cellEDepTotal_()
    , cellEDepIonizing_()
    , numCellPlane2d_()
    , numCellPlane3d_()
  {}

  //================================================================
  void MuCapSimHitsHist::beginJob() {

      art::ServiceHandle<art::TFileService> tfs;

      cellEDepTotal_ =  tfs->make<TH1D>("cellEDepTotal", "Cell energy deposition, total", 2000, 0., 2.);
      cellEDepTotal_->GetXaxis()->SetTitle("Energy deposition [MeV]");

      cellEDepIonizing_ =  tfs->make<TH1D>("cellEDepIonizing", "Cell energy deposition, ionizing", 2000, 0., 2.);
      cellEDepIonizing_->GetXaxis()->SetTitle("Energy deposition [MeV]");

      const int maxHitCells = 100;
      numCellPlane2d_ =  tfs->make<TH2D>("numCellPlane2d", "Num hit planes vs num hit cells",
                                         maxHitCells+1, -0.5, maxHitCells+0.5,
                                         57, -0.5, 56.5);

      numCellPlane2d_->SetOption("colz");


      numCellPlane3d_ =  tfs->make<TH3D>("numCellPlane3d", "Ek:Num hit planes:num hit cells",
                                         maxHitCells+1, -0.5, maxHitCells+0.5,
                                         29, -0.5, 28.5,
                                         pset_.get<unsigned>("numEnergyBins"),
                                         pset_.get<double>("ekMin"),
                                         pset_.get<double>("ekMax")
                                         );

      numCellPlane3d_->GetXaxis()->SetTitle("cells");
      numCellPlane3d_->GetYaxis()->SetTitle("planes");
      numCellPlane3d_->GetZaxis()->SetTitle("Ek, MeV");
  }

  //================================================================
  void MuCapSimHitsHist::analyze(const art::Event& event) {
    art::Handle<StepPointMCCollection> ih;
    event.getByLabel(hitsModuleLabel_, hitsInstanceName_, ih);
    const StepPointMCCollection& coll(*ih);
    //std::cout<<"MuCapSimHitsHist: got "<<coll.size()<<" sim hits"<<std::endl;

    typedef std::set<WirePlaneId> PlaneSet;
    PlaneSet hitPlanes;
    typedef std::map<WireCellId, double> CellMap;

    CellMap cellETotal;
    CellMap cellEIonizing;

    for(unsigned i=0; i<coll.size(); ++i) {
      const MuCapSimHit sh(coll[i]);
      // std::cout<<sh<<std::endl;

      hitPlanes.insert(sh.cid().plane());
      cellETotal[sh.cid()] += sh.hit().totalEDep();
      cellEIonizing[sh.cid()] += sh.hit().ionizingEdep();
    }

    const double primaryEk = getPrimaryEk(event);

    numCellPlane2d_->Fill(cellETotal.size(), hitPlanes.size());
    numCellPlane3d_->Fill(cellETotal.size(), hitPlanes.size(), primaryEk);

    for(CellMap::const_iterator i = cellETotal.begin(); i!=cellETotal.end(); ++i) {
      cellEDepTotal_->Fill(i->second);
    }
    for(CellMap::const_iterator i = cellEIonizing.begin(); i!=cellEIonizing.end(); ++i) {
      cellEDepIonizing_->Fill(i->second);
    }
  }

  //================================================================
  double MuCapSimHitsHist::getPrimaryEk(const art::Event& event) const {
    art::Handle<mu2e::SimParticleCollection> ih;
    event.getByLabel(particlesModuleLabel_, particlesInstanceName_, ih);

    cet::map_vector_key iprim(1);
    const mu2e::SimParticle primary(ih->getOrThrow(iprim));

    return primary.startMomentum().e() - primary.startMomentum().m();
  }

  //================================================================
} // namespace mucap

DEFINE_ART_MODULE(mucap::MuCapSimHitsHist);
