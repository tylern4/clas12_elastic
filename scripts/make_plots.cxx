#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <string>
#include "../include/constants.hpp"

int make_plots(std::string inFileName = "/Users/tylern/Desktop/show/clas12/test.root") {
  gStyle->SetHistMinimumZero();
  TFile *root_data = new TFile(inFileName.c_str());

  for (auto &&det : detector_name) {
    TCanvas *can = new TCanvas(Form("Phi_%s", det.second.c_str()), Form("Phi_%s", det.second.c_str()), 1600, 800);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *phi = (TH1D *)root_data->Get(Form("Phi/Phie_minus_Phip_%s_%d", det.second.c_str(), sec));
      phi->GetXaxis()->SetRangeUser(2, 4);
      phi->Draw("");
      TLine *l = new TLine(phi_min_cut, 0, phi_min_cut, phi->GetMaximum() * 1.05);
      l->Draw("SAME");
      TLine *r = new TLine(phi_max_cut, 0, phi_max_cut, phi->GetMaximum() * 1.05);
      r->Draw("SAME");
    }
  }

  {
    TCanvas *can = new TCanvas("MissingMass", "MissingMass", 1600, 800);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *missingMass2 = (TH1D *)root_data->Get(Form("W_vs_Q2/MM2_hist_sec_%d", sec));
      missingMass2->GetXaxis()->SetRangeUser(-1, 1);
      missingMass2->Draw("");
      TLine *l = new TLine(-MM2_cut, 0, -MM2_cut, missingMass2->GetMaximum() * 1.05);
      l->Draw("SAME");
      TLine *r = new TLine(MM2_cut, 0, MM2_cut, missingMass2->GetMaximum() * 1.05);
      r->Draw("SAME");
    }
  }

  for (auto &&det : detector_name) {
    TCanvas *can =
        new TCanvas(Form("W_at180_%s", det.second.c_str()), Form("W_at180_%s", det.second.c_str()), 1600, 800);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W = (TH1D *)root_data->Get(Form("at180/W_1pos_at180_%s_%d", det.second.c_str(), sec));
      W->Draw("");
    }
  }

  for (auto &&det : detector_name) {
    TCanvas *can =
        new TCanvas(Form("W_at180_MM_%s", det.second.c_str()), Form("W_at180_MM_%s", det.second.c_str()), 1600, 800);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W = (TH1D *)root_data->Get(Form("at180_MM/W_1pos_at180_MM_%s_%d", det.second.c_str(), sec));
      W->Draw("");
    }
  }

  for (auto &&det : detector_name) {
    TCanvas *can = new TCanvas(Form("W_at180_MMvsWo_%s", det.second.c_str()),
                               Form("W_at180_MMvsWo_%s", det.second.c_str()), 1600, 800);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W_wo = (TH1D *)root_data->Get(Form("at180/W_1pos_at180_%s_%d", det.second.c_str(), sec));
      W_wo->SetLineColor(kRed);
      W_wo->Draw("same");

      TH1D *W = (TH1D *)root_data->Get(Form("at180_MM/W_1pos_at180_MM_%s_%d", det.second.c_str(), sec));
      W->Draw("same");
    }
  }

  {
    TCanvas *can = new TCanvas("MomVsTheta", "MomVsTheta", 1600, 800);
    can->Divide(2, 2);
    can->cd(1);
    TH2D *MomVsTheta = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_both_0");
    MomVsTheta->Draw("");
    can->cd(4);
    TH2D *MomVsTheta_for = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_in_Forward_0");
    MomVsTheta_for->Draw("");
    can->cd(3);
    TH2D *MomVsTheta_cent = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_in_Central_0");
    MomVsTheta_cent->Draw("");
  }

  return 0;
}