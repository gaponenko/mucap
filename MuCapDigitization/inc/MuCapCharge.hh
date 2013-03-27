// Andrei Gaponenko, 2013

#ifndef MuCapDigitization_inc_MuCapCharge_hh
#define MuCapDigitization_inc_MuCapCharge_hh

#include <queue>
#include <map>

#include "MuCapDataProducts/inc/WireReadoutId.hh"

namespace mucap {

  struct MuCapTimedChargeDeposit {
    double time; // drift time is included
    double energy; // G4 deposit
    MuCapTimedChargeDeposit(double t, double e) : time(t), energy(e) {}

    // We accumulated deposits in a priority_queue, and want earlier times
    // to come out first.  Thus the inverted less-than definition:
    bool operator<(const MuCapTimedChargeDeposit& b) const {
      return b.time < this->time;
    }
  };

  // Queue ordered by time
  typedef std::priority_queue<MuCapTimedChargeDeposit> MuCapChargeHistory;

  typedef std::map<WireReadoutId,MuCapChargeHistory> MuCapChargeCollection;

}

#endif/*MuCapDigitization_inc_MuCapCharge_hh*/
