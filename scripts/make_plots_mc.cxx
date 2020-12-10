#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <THnSparse.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include "constants.hpp"

int make_plots_mc(std::string inFileName = "/Users/tylern/Desktop/show/clas12/test.root",
                  std::string mcFileName = "/Users/tylern/Desktop/show/clas12/test.root") {
  gStyle->SetHistMinimumZero();
  TFile *root_data = new TFile(inFileName.c_str());
  TFile *root_mc = new TFile(mcFileName.c_str());

  {
    TCanvas *can = new TCanvas("MomVsTheta_dataVsMC", "MomVsTheta_dataVsMC", 1920, 1080);
    can->Divide(2, 2);
    can->cd(1);
    TH2D *MomVsTheta = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_both_0");
    MomVsTheta->Draw("");
    can->cd(2);
    TH2D *MomVsTheta_mc = (TH2D *)root_mc->Get("Mom Vs Beta/MomVsTheta_pos_both_0");
    MomVsTheta_mc->Draw("");
    can->cd(3);
    TH2D *MomVsTheta_cent = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_in_Central_0");
    MomVsTheta_cent->Draw("");
    can->cd(4);
    TH2D *MomVsTheta_for = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_in_Forward_0");
    MomVsTheta_for->Draw("");
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  {
    TCanvas *can = new TCanvas("MomVsTheta_cutLow", "MomVsTheta_cutLow", 1920, 1080);
    can->Divide(2, 2);
    can->cd(1);
    TH2D *MomVsTheta = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_both_0");
    MomVsTheta->Draw("");
    can->cd(2);
    TH2D *MomVsTheta_mc = (TH2D *)root_mc->Get("Mom Vs Beta/MomVsTheta_pos_both_0");
    MomVsTheta_mc->Draw("");
    can->cd(3);
    TH2D *MomVsTheta_cent = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_in_Central_0");
    MomVsTheta_cent->Draw("");
    can->cd(4);
    TH2D *MomVsTheta_for = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_in_Forward_0");
    MomVsTheta_for->Draw("");
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  return 0;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = make_plots_mc(argv[1]);

  return can;
}
#endif