// Andrei Gaponenko, 2012

#ifndef postproc_mecofitter_applyHypothesis_hh
#define postproc_mecofitter_applyHypothesis_hh

#include <cassert>

#include "FitData.hh"
#include "FitMC.hh"

namespace mucap {
  namespace fitter {

    // SimpleHypothesis is a proton spectrum shape with no
    // free parameters (e.g. a double=>double function).
    template<class SimpleHypothesis>
    FitData applyHypothesis(const FitMC& mc, const SimpleHypothesis& spectrum) {
      assert(!mc.empty());

      FitData res(mc.front().values.size(), 0.);

      for(FitMC::const_iterator iz = mc.begin(); iz != mc.end(); ++iz) {
        const double w = spectrum(iz->energy);
        for(unsigned i=0; i<res.size(); ++i) {
          res[i] += w * iz->values[i];
        }
      }

      return res;
    }

  } // fitter
} // mucap

#endif/*postproc_mecofitter_applyHypothesis_hh*/
