// Andrei Gaponenko, 2013

#include "MuCapGenerator/inc/SpectrumGenTabulated.hh"
#include "MuCapUtilities/inc/TabulatedFunction.hh"

#include <cassert>

namespace mucap {

  std::vector<double> SpectrumGenTabulated::getYValues(const TabulatedFunction& func) {
    std::vector<double> res;
    for(const auto& p : func.table()) {
      res.emplace_back(p.y);
    }
    return res;
  }

  SpectrumGenTabulated::SpectrumGenTabulated(const TabulatedFunction& func, CLHEP::HepRandomEngine& rng)
    : tmp_(getYValues(func))
    , spectrum_(rng, &tmp_[0], tmp_.size())
    , xmin_(func.table().front().x)
    , xmax_(func.table().back().x)
  {
    assert(func.table().size() > 1);
    tmp_.clear();
  }

  double SpectrumGenTabulated::generate() {
    return xmin_ + (xmax_ - xmin_)*spectrum_.fire();
  }

} // namespace mucap
