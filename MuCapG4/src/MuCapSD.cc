// Andrei Gaponenko, 2012

#include "MuCapG4/inc/MuCapSD.hh"

#include <sstream>

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

namespace mucap {

  using namespace mu2e;

  //================================================================
  const std::string& MuCapSD::name() {
    static std::string name_("MuCapSimHits");
    return name_;
  };

  //================================================================
  MuCapSD::MuCapSD(const fhicl::ParameterSet& pset)
    : G4VSensitiveDetector(name())
    , processInfo_()
    , collection_()
    , simID_()
    , event_()
    , sizeLimit_(pset.get<unsigned>("hitsSizeLimit"))
    , totalHitCount_(0)
  {}

  //================================================================
  void MuCapSD::beforeG4Event(StepPointMCCollection *outputHits,
                              PhysicsProcessInfo *processInfo,
                              const art::ProductID& simID,
                              const art::Event& event )
  {
    processInfo_ = processInfo;
    collection_  = outputHits;
    simID_       = &simID;
    event_       = &event;
    totalHitCount_ = 0;
  }

  //================================================================
  G4bool MuCapSD::ProcessHits(G4Step* aStep,G4TouchableHistory*) {

    bool retval = false;

    G4double totalEDep = aStep->GetTotalEnergyDeposit();

    if(totalEDep > 0.) {

      ++totalHitCount_;

      if(!sizeLimit_||(collection_->size() < sizeLimit_)) {

        retval = true;

        art::Ptr<SimParticle> particle( *simID_, aStep->GetTrack()->GetTrackID(), event_->productGetter(*simID_) );

        // Which process caused this step to end?
        G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        ProcessCode endCode(processInfo_->findAndCount(pname));

        collection_->
          push_back(StepPointMC(particle,
                                aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                                aStep->GetTotalEnergyDeposit(),
                                aStep->GetNonIonizingEnergyDeposit(),
                                aStep->GetPreStepPoint()->GetGlobalTime(),
                                aStep->GetPreStepPoint()->GetProperTime(),
                                aStep->GetPreStepPoint()->GetPosition(),
                                aStep->GetPreStepPoint()->GetMomentum(),
                                aStep->GetStepLength(),
                                endCode
                                ));
      }
    }

    return retval;
  }

  //================================================================
  void MuCapSD::EndOfEvent(G4HCofThisEvent*){

    if(sizeLimit_>0 && (collection_->size() >= sizeLimit_)){
      std::ostringstream os;
      os<< "WARNING: " << totalHitCount_<< " MuCap hits were generated in the event.\n"
        << "Only " <<collection_->size()<< " are saved in output collection.\n";

      mf::LogWarning("G4")<<os.str();
      G4cout<<os.str()<<G4endl;
    }

  }

  //================================================================

} //namespace mucap
