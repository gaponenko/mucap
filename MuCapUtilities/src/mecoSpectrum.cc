#include "MuCapUtilities/inc/mecoSpectrum.hh"
#include <cmath>

double mucap::mecoSpectrum(double T, const MECOPars& pars) {
  return (T > pars.Tth) ?
    pars.A*std::pow((1-pars.Tth/T), pars.alpha) * exp(-T/pars.T0) :
    0.;
}

//================================================================
// g++  -o mecoSpectrum -DBUILD_TEST -I. MuCapUtilities/src/mecoSpectrum.cc
#ifdef BUILD_TEST
#include <iostream>
int main() {
  using namespace mucap;
  MECOPars p;
  p.A = 0.105;
  p.Tth = 1.4;
  p.alpha = 1.3279;
  p.T0 = 3.1;

  for(double ek = 0; ek < 100; ek += 0.5) {
    std::cout<<ek<<"\t"<<mecoSpectrum(ek, p)<<std::endl;
  }

  return 0;
}

#endif
//================================================================
