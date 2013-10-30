// Andrei Gaponenko, 2012

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iomanip>

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
  class MuCapG3CmpPrinter : public art::EDAnalyzer {
    std::string inModuleLabel_;
    std::string inInstanceName_;
    std::string outFileName_;
    std::ofstream of_;
  public:
    explicit MuCapG3CmpPrinter(const fhicl::ParameterSet& pset);
    virtual void beginJob() override;
    virtual void analyze(const art::Event& event) override;
  };

  //================================================================
  MuCapG3CmpPrinter::MuCapG3CmpPrinter(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
    , inInstanceName_(pset.get<std::string>("particlesInstanceName", ""))
    , outFileName_(pset.get<std::string>("outFileName"))
  {}

  //================================================================
  void MuCapG3CmpPrinter::beginJob() {
    //art::ServiceHandle<art::TFileService> tfs;
    of_.open(outFileName_.c_str());
    if(!of_.is_open()) {
      throw cet::exception("RUNTIME")<<"Could not open file "<<outFileName_<<" for writing\n";
    }
  }

  //================================================================
  void MuCapG3CmpPrinter::analyze(const art::Event& event) {

    art::Handle<SimParticleCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
    const SimParticleCollection& coll(*ih);
    
    cet::map_vector_key iprim(1);
    const SimParticle primary(coll.getOrThrow(iprim));
    if(primary.pdgId() != 2212) {
      throw cet::exception("INPUTS")<<"MuCapG3CmpPrinter: assumed primary is not a proton: pdgId = "<<primary.pdgId()<<"\n";
    }

    //std::cout<<primary.startPosition()<<" sm="<<primary.startMomentum()<<", ep="<<primary.endPosition()<< ", em="<<primary.endMomentum()<<std::endl;

    // Immitate the fortran printout from G3
    of_<<std::scientific<<std::setprecision(8)
       <<std::setw(15)<<primary.startPosition().x()<<" "
       <<std::setw(15)<<primary.startPosition().y()<<" "
       <<std::setw(15)<<primary.startPosition().z()<<" "
       <<std::setw(15)<<primary.startMomentum().x()<<" "
       <<std::setw(15)<<primary.startMomentum().y()<<" "
       <<std::setw(15)<<primary.startMomentum().z()<<" "
       <<std::setw(15)<<primary.endPosition().x()<<" "
       <<std::setw(15)<<primary.endPosition().y()<<" "
       <<std::setw(15)<<primary.endPosition().z()
       <<std::endl;
  }

  //================================================================
} // namespace mucap

DEFINE_ART_MODULE(mucap::MuCapG3CmpPrinter);
