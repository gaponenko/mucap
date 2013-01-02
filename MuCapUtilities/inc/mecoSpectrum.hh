// An empirical function used to fit to the ejected proton spectrum in
// MECO note 034 (Mu2e docdb-105) by Ed V. Hungerford.
// The fitted parameter values are:
//
// A = 0.105/MeV, Tth=1.4 MeV, alpha=1.3279, T0=3.1 MeV
//
// We use T0inv == 1/T0 instead of the T0 parameter as this is
// more convenient for making the fit.
//
// Andrei Gaponenko, 2012

#ifndef MuCapUtilities_inc_mecoSpectrum_hh
#define MuCapUtilities_inc_mecoSpectrum_hh

#include <ostream>

namespace mucap {
  struct MECOPars {
    double A; // normalization
    double Tth; // threshold
    double alpha;
    double T0inv;
    MECOPars() : A(), Tth(), alpha(), T0inv() {}
  };

  std::ostream& operator<<(std::ostream& os, const MECOPars& p);

  double mecoSpectrum(double ek, const MECOPars& pars);

}
#endif/*MuCapUtilities_inc_mecoSpectrum_hh*/
