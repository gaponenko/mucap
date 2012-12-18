// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class MuCapSimHitsHist : public art::EDAnalyzer {
    std::string inModuleLabel_;
    std::string inInstanceName_;

    //EMFSimHitHistograms ch_;

  public:
    explicit MuCapSimHitsHist(const fhicl::ParameterSet& pset);
    //virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  MuCapSimHitsHist::MuCapSimHitsHist(const fhicl::ParameterSet& pset)
    : inModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inInstanceName_(pset.get<std::string>("inputInstanceName", ""))
  {}

  //================================================================
  void MuCapSimHitsHist::analyze(const art::Event& event) {
    art::Handle<StepPointMCCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
    const StepPointMCCollection& coll(*ih);
    std::cout<<"MuCapSimHitsHist: got "<<coll.size()<<" sim hits"<<std::endl;
    for(unsigned i=0; i<coll.size(); ++i) {
      const StepPointMC& hit = coll[i];
      std::cout<<"Hit vid="<<hit.volumeId()<<", edep = "<<hit.totalEDep()<<std::endl;
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuCapSimHitsHist);
