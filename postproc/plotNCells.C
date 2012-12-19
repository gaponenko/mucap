#include "getHisto.h"

TH1 *getNCells(std::string file) {
  TH2 *h2 = dynamic_cast<TH2*>(getHisto(file, "mchits/numCellPlane2d"));
  TH1 *h1 = h2->ProjectionX();
  h1->SetDirectory(0);

  h1->GetXaxis()->SetRangeUser(-0.5, 49.5);
  h1->GetXaxis()->SetTitle("Number of hit wires");

  return h1;
}

void plotNCells() {
  TCanvas *cc = makeEpsCanvas("numHitWires");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH1 *h25 = getNCells("mcanalysis_2.5.root");
  TH1 *h50 = getNCells("mcanalysis_5.0.root");
  TH1 *h75 = getNCells("mcanalysis_7.5.root");

  h25->SetMarkerStyle(22);
  h25->SetMarkerColor(kRed);
  h25->SetLineColor(kRed);

  h50->SetMarkerStyle(21);
  h50->SetMarkerColor(kMagenta);
  h50->SetLineColor(kMagenta);

  h75->SetMarkerStyle(20);
  h75->SetMarkerColor(kBlue);
  h75->SetLineColor(kBlue);

  h25->Draw("PE");
  h50->Draw("PE same");
  h75->Draw("PE same");

  TLegend *leg = new TLegend(0.68, 0.68, 0.88, 0.88);
  leg->SetName("legend");
  leg->AddEntry(h75, "7.5 MeV", "p");
  leg->AddEntry(h50, "5.0 MeV", "p");
  leg->AddEntry(h25, "2.5 MeV", "p");
  leg->Draw();

  cc->Print(".pdf");
}
