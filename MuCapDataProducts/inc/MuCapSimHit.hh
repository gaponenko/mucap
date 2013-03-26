// This class does not have to be persisted because it can be easily re-created from StepPointMC.
//
// Andrei Gaponenko, 2012

#ifndef MuCapDataProducts_inc_MuCapSimHit_hh
#define MuCapDataProducts_inc_MuCapSimHit_hh

#include <ostream>

#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mu2e { class StepPointMC; }

namespace mucap {

  class MuCapSimHit {
  public:

    explicit MuCapSimHit(const mu2e::StepPointMC& sp);

    const mu2e::StepPointMC& hit() const { return *sp_; }
    const WireCellId&  cid() const { return cid_; }

  private:
    const mu2e::StepPointMC *sp_;
    WireCellId cid_;

  };

  std::ostream& operator<<(std::ostream& os, const MuCapSimHit& sh);

}

#endif /* MuCapDataProducts_inc_MuCapSimHit_hh */
