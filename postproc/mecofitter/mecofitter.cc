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
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "Math/Minimizer.h"

#include "TemplateFitter.hh"

int main(int argc, const char* argv[]) {
  using fhicl::ParameterSet;
  using std::string;

  if(argc != 2) {
    std::cerr<<"Usage: mecofitter config.flc\n";
    exit(1);
  }

  int ret = 1;

  try {
    ParameterSet pset;
    cet::filepath_maker pm;
    fhicl::make_ParameterSet(argv[1], pm, pset);

    TH1::AddDirectory(kFALSE);

    //----------------------------------------------------------------
    mucap::fitter::TemplateFitter tf(pset);

    std::cout<<"Minimization status = "<<tf.getMinimizer().Minimize()<<std::endl;

    const double* px = tf.getMinimizer().X();
    std::cout<<"Parameter values = "
             <<px[0]<<", "
             <<px[1]<<", "
             <<px[2]<<", "
             <<px[3]
             <<std::endl;

    std::cout<<std::endl;
    tf.getMinimizer().PrintResults();
    std::cout<<std::endl;

    //----------------------------------------------------------------
    // Write stuff to the output ROOT file

    const string outfilename(pset.get<ParameterSet>("outputs").get<string>("filename"));
    TFile out(outfilename.c_str(), "RECREATE");
    if(!out.IsOpen()) {
      throw cet::exception("INPUTS")
        <<"Error opening output file "<<outfilename<<"\n";
    }

    tf.visualizeFitRegion()->Write();
    tf.visualizeFitData("data", tf.getFunction().getFitData())->Write();

    out.Close();

    ret = 0;
  }
  catch(std::exception& e) {
    std::cerr<<"Got std::exception: "<<e.what()<<std::endl;
  }

  return ret;
}
