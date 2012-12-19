#include "getHisto.h"

void plot2d(std::string file) {
  TH2 *h2 = dynamic_cast<TH2*>(getHisto(file, "mchits/numCellPlane2d"));
  h2->GetXaxis()->SetRangeUser(-0.5, 49.5);
  h2->GetXaxis()->SetTitle("Number of hit wires");

  h2->GetYaxis()->SetRangeUser(-0.5, 35.5);
  h2->GetYaxis()->SetTitle("Number of hit planes");

  TCanvas *cc = makeEpsCanvas("planesVsWires");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  h2->Draw();
}
