// Andrei Gaponenko, 2012

#ifndef MuCapDataProducts_inc_WirePlaneId_hh
#define MuCapDataProducts_inc_WirePlaneId_hh

#include <ostream>

namespace mu2e {

  class WirePlaneId {
  public:

    static const unsigned int NOPLANE = -1u;

    // No automatic conversion of int to WirePlaneId.
    explicit WirePlaneId(unsigned int plane) : plane_(plane) {}

    // zero based
    unsigned int number() const { return plane_;}

    bool operator==( WirePlaneId const& rhs) const{
      return (plane_ == rhs.plane_);
    }

    bool operator!=( WirePlaneId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( WirePlaneId const& rhs) const{
      return (plane_ < rhs.plane_);
    }

    // Default constructor is required by ROOT persistency
    WirePlaneId() : plane_(NOPLANE) {}

  private:
    unsigned int plane_;
  };

  inline std::ostream& operator<<( std::ostream& os, const WirePlaneId& id) {
    return os<<id.number();
  }
}
#endif /* MuCapDataProducts_inc_WirePlaneId_hh */
