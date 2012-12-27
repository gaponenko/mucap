// Fitting a hit multiplicity distribution by varying parameters of
// the spectrum of ejected protons.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "cetlib/filepath_maker.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "TH1.h"
#include "TH3.h"

#include "TemplateFitter.hh"

int main(int argc, const char* argv[]) {
  if(argc != 2) {
    std::cerr<<"Usage: mecofitter config.flc\n";
    exit(1);
  }

  int ret = 1;

  try {
    //cet::filepath_lookup pm("FHICL_FILE_PATH");
    cet::filepath_maker pm;
    fhicl::ParameterSet pset;
    fhicl::make_ParameterSet(argv[1], pm, pset);

    TH1::AddDirectory(kFALSE);

    mucap::fitter::TemplateFitter tf(pset);

    ret = 0;
  }
  catch(std::exception& e) {
    std::cerr<<"Got std::exception: "<<e.what()<<std::endl;
  }

  return ret;
}
