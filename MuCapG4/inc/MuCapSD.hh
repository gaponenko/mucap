// Andrei Gaponenko, 2012

#ifndef MuCapG4_MuCapSD_hh
#define MuCapG4_MuCapSD_hh

#include "fhiclcpp/ParameterSet.h"

#include "art/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

#include "G4VSensitiveDetector.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  class PhysicsProcessInfo;

  // This class should not derive from Mu2eSensitiveDetector because the latter
  // presumes a wrong type for the hit collection.
  class MuCapSD : public G4VSensitiveDetector {

    // Non-ownning pointer and object that returns code describing physics processes.
    PhysicsProcessInfo* processInfo_;

    // Non-owning pointer to the  collection into which hits will be added.
    StepPointMCCollection* collection_;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    const art::ProductID *simID_;
    const art::Event *event_;

    // Limit maximum size of the steps collection
    unsigned sizeLimit_;
    unsigned totalHitCount_;

  public:

    static const std::string& name();

    explicit MuCapSD(const fhicl::ParameterSet& pset);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection *outputHits,
                       PhysicsProcessInfo *processInfo,
                       const art::ProductID& simID,
                       const art::Event& event);
  };

} // namespace mu2e

#endif /* MuCapG4_MuCapSD_hh */
