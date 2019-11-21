/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");
  q2_max = 10.6;
  w_max = 4.5;
  p_max = 10.6;

  makeHists();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  Write_WvsQ2();

  std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();

  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  short pos_det = _e->pos_det();
  if ((sec > 0 && sec < num_sectors) || pos_det != -1) {
    W_hist_all_events[0]->Fill(_e->W());
    W_vs_q2_all_events[0]->Fill(_e->W(), _e->Q2());
    W_hist_all_events[sec]->Fill(_e->W());
    W_vs_q2_all_events[sec]->Fill(_e->W(), _e->Q2());

    if (_e->onePositive()) {
      W_hist_1pos[0]->Fill(_e->W());
      W_vs_q2_1pos[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos[sec]->Fill(_e->W());
      W_vs_q2_1pos[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_MM0()) {
      W_hist_1pos_0charge[0]->Fill(_e->W());
      W_vs_q2_1pos_0charge[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_0charge[sec]->Fill(_e->W());
      W_vs_q2_1pos_0charge[sec]->Fill(_e->W(), _e->Q2());

      Phie_vs_Phip[pos_det][0]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][0]->Fill(_e->phi_diff());
      Phie_vs_Phip[pos_det][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][sec]->Fill(_e->phi_diff());
    }
    if (_e->onePositive_part2()) {
      W_hist_1pos_gpart2[0]->Fill(_e->W());
      W_vs_q2_1pos_gpart2[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_gpart2[sec]->Fill(_e->W());
      W_vs_q2_1pos_gpart2[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at90()) {
      MissingMass[0]->Fill(_e->MM2());
      MissingMass[sec]->Fill(_e->MM2());
    }
    if (_e->onePositive_at90_MM0()) {
      W_hist_1pos_at90[pos_det][0]->Fill(_e->W());
      W_vs_q2_1pos_at90[pos_det][0]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at90[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at90[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
  }
}

void Histogram::Write_WvsQ2() {
  TDirectory* phi_folder = RootOutputFile->mkdir("Phi");
  phi_folder->cd();
  for (short j = 0; j < 2; j++) {
    for (int i = 0; i < num_sectors; i++) {
      Phie_vs_Phip[j][i]->SetXTitle("Phie");
      Phie_vs_Phip[j][i]->SetYTitle("Phip");
      Phie_vs_Phip[j][i]->SetOption("COLZ");
      Phie_vs_Phip[j][i]->Write();
    }
    for (int i = 0; i < num_sectors; i++) {
      Phie_Phip_hist[j][i]->SetXTitle("Phi");
      Phie_Phip_hist[j][i]->Write();
    }
  }

  TDirectory* at90_folder = RootOutputFile->mkdir("at90");
  at90_folder->cd();
  for (short j = 0; j < 2; j++) {
    for (int i = 0; i < num_sectors; i++) {
      W_hist_1pos_at90[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at90[j][i]->Write();
    }
    for (int i = 0; i < num_sectors; i++) {
      W_vs_q2_1pos_at90[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at90[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at90[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at90[j][i]->Write();
    }
  }

  TDirectory* W_vs_Q2_folder = RootOutputFile->mkdir("W_vs_Q2");
  W_vs_Q2_folder->cd();
  for (int i = 0; i < num_sectors; i++) {
    MissingMass[i]->SetXTitle("MM^2 (GeV)");
    MissingMass[i]->Write();

    W_hist_all_events[i]->SetXTitle("W (GeV)");
    W_hist_all_events[i]->Write();
    W_hist_1pos[i]->SetXTitle("W (GeV)");
    W_hist_1pos[i]->Write();
    W_hist_1pos_0charge[i]->SetXTitle("W (GeV)");
    W_hist_1pos_0charge[i]->Write();
    W_hist_1pos_gpart2[i]->SetXTitle("W (GeV)");
    W_hist_1pos_gpart2[i]->Write();

    W_vs_q2_all_events[i]->SetXTitle("W (GeV)");
    W_vs_q2_all_events[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_all_events[i]->SetOption("COLZ");
    W_vs_q2_all_events[i]->Write();

    W_vs_q2_1pos[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos[i]->SetOption("COLZ");
    W_vs_q2_1pos[i]->Write();

    W_vs_q2_1pos_0charge[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_0charge[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_0charge[i]->SetOption("COLZ");
    W_vs_q2_1pos_0charge[i]->Write();

    W_vs_q2_1pos_gpart2[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_gpart2[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_gpart2[i]->SetOption("COLZ");
    W_vs_q2_1pos_gpart2[i]->Write();
  }
}

void Histogram::makeHists() {
  for (short sec = 0; sec < num_sectors; sec++) {
    MissingMass[sec] =
        std::make_shared<TH1D>(Form("MM2_hist_sec_%d", sec), Form("MM2_hist_sec_%d", sec), bins, -w_max, w_max);
    W_hist_all_events[sec] =
        std::make_shared<TH1D>(Form("W_hist_sec_%d", sec), Form("W_hist_sec_%d", sec), bins, zero, w_max);
    W_hist_1pos[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_%d", sec), Form("W_hist_1pos_%d", sec), bins, zero, w_max);
    W_hist_1pos_0charge[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_MM0_%d", sec), Form("W_hist_1pos_MM0_%d", sec), bins, zero, w_max);
    W_hist_1pos_gpart2[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_part2_%d", sec), Form("W_hist_1pos_part2_%d", sec), bins, zero, w_max);

    W_vs_q2_all_events[sec] =
        std::make_shared<TH2D>(Form("WQ2_sec_%d", sec), Form("WQ2_sec_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_%d", sec), Form("WQ2_1pos_%d", sec), bins, zero, w_max,
                                               bins, zero, q2_max);
    W_vs_q2_1pos_0charge[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_MM0_%d", sec), Form("WQ2_1pos_MM0_%d", sec), bins,
                                                       zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos_gpart2[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_part2_%d", sec), Form("WQ2_1pos_part2_%d", sec),
                                                      bins, zero, w_max, bins, zero, q2_max);
    for (auto&& det : detector_name) {
      int d = 0;
      if (det.first == 2) d = 0;
      if (det.first == 4) d = 1;
      MomVsBeta[d][sec] =
          std::make_shared<TH2D>(Form("MomVsBeta_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsBeta_%s_%d", det.second.c_str(), sec), 500, zero, p_max, 500, zero, 1.2);
      Phie_vs_Phip[d][sec] =
          std::make_shared<TH2D>(Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec), 500, -PI, PI, 500, -PI, PI);
      Phie_Phip_hist[d][sec] =
          std::make_shared<TH1D>(Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec), 500, zero, 2 * PI);
      W_hist_1pos_at90[d][sec] =
          std::make_shared<TH1D>(Form("W_1pos_at90_%s_%d", det.second.c_str(), sec),
                                 Form("W_1pos_at90_%s_%d", det.second.c_str(), sec), bins, zero, w_max);
      W_vs_q2_1pos_at90[d][sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at90_%s_%d", det.second.c_str(), sec),
                                                         Form("WQ2_1pos_at90_%s_%d", det.second.c_str(), sec), bins,
                                                         zero, w_max, bins, zero, q2_max);
    }
  }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e) {
  MomVsBeta[_e->pos_det()][0]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());
}

void Histogram::Write_MomVsBeta() {
  for (short i = 0; i < 2; i++)
    for (short p = 0; p < num_sectors; p++) {
      MomVsBeta[i][p]->SetXTitle("Momentum (GeV)");
      MomVsBeta[i][p]->SetYTitle("#beta");
      MomVsBeta[i][p]->SetOption("COLZ1");
      MomVsBeta[i][p]->Write();
    }
}
