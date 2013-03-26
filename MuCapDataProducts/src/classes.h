#include "MuCapDataProducts/inc/WireReadoutId.hh"
#include "MuCapDataProducts/inc/MuCapRecoHit.hh"
#include "MuCapDataProducts/inc/MuCapRecoHitCollection.hh"

#include "art/Persistency/Common/Wrapper.h"

#include <vector>

template class std::vector<mucap::MuCapRecoHit>;
template class art::Wrapper<mucap::MuCapRecoHitCollection>;
