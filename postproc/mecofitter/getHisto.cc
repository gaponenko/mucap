#include "getHisto.hh"

#include <stdexcept>
#include "TFile.h"
#include "TH1.h"

//================================================================
TH1* getHisto(const std::string& filename, const std::string& histoname) {
  TFile fff(filename.c_str());
  if(!fff.IsOpen()) {
    throw std::runtime_error("getHisto(): Error opening file \""+filename+"\"");
  }
  TH1* h = dynamic_cast<TH1*>(fff.Get(histoname.c_str()));
  if(!h) {
    throw std::runtime_error("No such histogram \""+histoname
                             +"\" in file \""+filename+"\"");
  }
  h->SetDirectory(0);
  return h;
}

//================================================================
