// Andrei Gaponenko, 2012

#ifndef MuCapDataProducts_inc_MuCapRecoHit_hh
#define MuCapDataProducts_inc_MuCapRecoHit_hh

#include <ostream>

#include "MuCapDataProducts/inc/WireReadoutId.hh"

namespace mucap {

  class MuCapRecoHit {
  public:

    MuCapRecoHit() : time_(), width_() {}
    MuCapRecoHit(const WireReadoutId& rid, double t, double w) : rid_(rid), time_(t), width_(w) {}

    const WireReadoutId&  rid() const { return rid_; }
    double time() const { return time_; }
    double width() const { return width_; }

  private:
    WireReadoutId rid_;
    double time_;
    double width_;
  };

  std::ostream& operator<<(std::ostream& os, const MuCapRecoHit& sh);

}

#endif /* MuCapDataProducts_inc_MuCapRecoHit_hh */
