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
#include "MuCapDataProducts/inc/WirePlaneId.hh"
#include "MuCapDataProducts/inc/WireCellId.hh"
#include "MuCapDataProducts/inc/MuCapSimHit.hh"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class MuCapSimHitsHist : public art::EDAnalyzer {
    std::string inModuleLabel_;
    std::string inInstanceName_;

    TH1 *cellEDepTotal_;
    TH1 *cellEDepIonizing_;
    TH2 *numCellPlane2d_;

  public:
    explicit MuCapSimHitsHist(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  MuCapSimHitsHist::MuCapSimHitsHist(const fhicl::ParameterSet& pset)
    : inModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inInstanceName_(pset.get<std::string>("inputInstanceName", ""))
    , cellEDepTotal_()
    , cellEDepIonizing_()
    , numCellPlane2d_()
  {}

  //================================================================
  void MuCapSimHitsHist::beginJob() {

      art::ServiceHandle<art::TFileService> tfs;

      cellEDepTotal_ =  tfs->make<TH1D>("cellEDepTotal", "Cell energy deposition, total", 2000, 0., 2.);
      cellEDepTotal_->GetXaxis()->SetTitle("Energy deposition [MeV]");

      cellEDepIonizing_ =  tfs->make<TH1D>("cellEDepIonizing", "Cell energy deposition, ionizing", 2000, 0., 2.);
      cellEDepIonizing_->GetXaxis()->SetTitle("Energy deposition [MeV]");

      const int maxHitCells = 300;
      numCellPlane2d_ =  tfs->make<TH2D>("numCellPlane2d", "Num hit planes vs num hit cells",
                                         maxHitCells+1, -0.5, maxHitCells+0.5,
                                         57, -0.5, 56.5);

      numCellPlane2d_->SetOption("colz");
  }

  //================================================================
  void MuCapSimHitsHist::analyze(const art::Event& event) {
    art::Handle<StepPointMCCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
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

    numCellPlane2d_->Fill(cellETotal.size(), hitPlanes.size());
    for(CellMap::const_iterator i = cellETotal.begin(); i!=cellETotal.end(); ++i) {
      cellEDepTotal_->Fill(i->second);
    }
    for(CellMap::const_iterator i = cellEIonizing.begin(); i!=cellEIonizing.end(); ++i) {
      cellEDepIonizing_->Fill(i->second);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuCapSimHitsHist);
