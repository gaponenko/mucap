// Andrei Gaponenko, 2012

#ifndef MuCapDataProducts_inc_WireCellId_hh
#define MuCapDataProducts_inc_WireCellId_hh

#include <ostream>

#include "MuCapDataProducts/inc/WirePlaneId.hh"

namespace mu2e {

  class WireCellId {
  public:

    WireCellId(const WirePlaneId& plane, unsigned int cell);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    WireCellId() : plane_(-1), cell_() {}

    const WirePlaneId& plane() const { return plane_; }
    unsigned int cell() const { return cell_; }

    unsigned encodeToInteger() const;
    static WireCellId decodeFromInteger(unsigned encodedVolumeId);

    bool operator==(const WireCellId& rhs) const {
      return (plane_ == rhs.plane_)&&(cell_ == rhs.cell_);
    }

    bool operator!=( WireCellId const& rhs) const {
      return !(*this == rhs);
    }

    bool operator<(const WireCellId& rhs) const {
      return
        (plane_ < rhs.plane_) ||
        ((plane_ == rhs.plane_) && ((cell_ < rhs.cell_) )
         );
    }

  private:
    WirePlaneId plane_;
    unsigned int cell_;
  };

  //----------------------------------------------------------------
  std::ostream& operator<<( std::ostream& os, const WireCellId& id);
}
#endif /* MuCapDataProducts_inc_WireCellId_hh */
