// Readouts correspond to drift volumes assuming the wires.
// Several readouts may be read out on the same channel - see WireReadoutId.
//
// Andrei Gaponenko, 2012

#ifndef MuCapDataProducts_inc_WireReadoutId_hh
#define MuCapDataProducts_inc_WireReadoutId_hh

#include <ostream>

#include "MuCapDataProducts/inc/WireCellId.hh"

namespace mucap {

  class WireReadoutId {
  public:

    WireReadoutId(const WirePlaneId& plane, unsigned int channel);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    WireReadoutId() : plane_(-1), channel_() {}

    const WirePlaneId& plane() const { return plane_; }
    unsigned int channel() const { return channel_; }

    bool operator==(const WireReadoutId& rhs) const {
      return (plane_ == rhs.plane_)&&(channel_ == rhs.channel_);
    }

    bool operator!=( WireReadoutId const& rhs) const {
      return !(*this == rhs);
    }

    bool operator<(const WireReadoutId& rhs) const {
      return
        (plane_ < rhs.plane_) ||
        ((plane_ == rhs.plane_) && ((channel_ < rhs.channel_) )
         );
    }

  private:
    WirePlaneId plane_;
    unsigned int channel_;
  };

  //----------------------------------------------------------------
  std::ostream& operator<<( std::ostream& os, const WireReadoutId& id);
}

#endif /* MuCapDataProducts_inc_WireReadoutId_hh */
