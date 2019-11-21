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
  if (sec != 0 && sec < num_sectors) {
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
    if (_e->onePositive_0Q()) {
      W_hist_1pos_0charge[0]->Fill(_e->W());
      W_vs_q2_1pos_0charge[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_0charge[sec]->Fill(_e->W());
      W_vs_q2_1pos_0charge[sec]->Fill(_e->W(), _e->Q2());

      Phie_vs_Phip[0]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[0]->Fill(_e->phi_diff());
      Phie_vs_Phip[sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[sec]->Fill(_e->phi_diff());
    }
    if (_e->onePositive_part2()) {
      W_hist_1pos_gpart2[0]->Fill(_e->W());
      W_vs_q2_1pos_gpart2[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_gpart2[sec]->Fill(_e->W());
      W_vs_q2_1pos_gpart2[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at90()) {
      W_hist_1pos_at90[0]->Fill(_e->W());
      W_vs_q2_1pos_at90[0]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at90[sec]->Fill(_e->W());
      W_vs_q2_1pos_at90[sec]->Fill(_e->W(), _e->Q2());
    }

  } else {
    std::cerr << "Sector " << sec << std::endl;
  }
}

void Histogram::Write_WvsQ2() {
  TDirectory* phi_folder = RootOutputFile->mkdir("Phi");
  phi_folder->cd();
  for (int i = 0; i < num_sectors; i++) {
    Phie_vs_Phip[i]->SetXTitle("Phie");
    Phie_vs_Phip[i]->SetYTitle("Phip");
    Phie_vs_Phip[i]->SetOption("COLZ");
    Phie_vs_Phip[i]->Write();
  }
  for (int i = 0; i < num_sectors; i++) {
    Phie_Phip_hist[i]->SetXTitle("Phi");
    Phie_Phip_hist[i]->Write();
  }

  TDirectory* at90_folder = RootOutputFile->mkdir("at90");
  at90_folder->cd();
  for (int i = 0; i < num_sectors; i++) {
    W_hist_1pos_at90[i]->SetXTitle("W (GeV)");
    W_hist_1pos_at90[i]->Write();
  }
  for (int i = 0; i < num_sectors; i++) {
    W_vs_q2_1pos_at90[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_at90[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_at90[i]->SetOption("COLZ");
    W_vs_q2_1pos_at90[i]->Write();
  }

  TDirectory* W_vs_Q2_folder = RootOutputFile->mkdir("W_vs_Q2");
  W_vs_Q2_folder->cd();
  for (int i = 0; i < num_sectors; i++) {
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
    MomVsBeta[sec] =
        std::make_shared<TH2D>(Form("MomVsBeta_%d", sec), Form("MomVsBeta_%d", sec), 500, zero, p_max, 500, zero, 1.2);
    Phie_vs_Phip[sec] =
        std::make_shared<TH2D>(Form("Phie_vs_Phip_%d", sec), Form("Phie_vs_Phip_%d", sec), 500, -PI, PI, 500, -PI, PI);
    Phie_Phip_hist[sec] =
        std::make_shared<TH1D>(Form("Phie_minus_Phip_%d", sec), Form("Phie_minus_Phip_%d", sec), 500, zero, 2 * PI);

    W_hist_all_events[sec] =
        std::make_shared<TH1D>(Form("W_hist_sec_%d", sec), Form("W_hist_sec_%d", sec), bins, zero, w_max);
    W_hist_1pos[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_%d", sec), Form("W_hist_1pos_%d", sec), bins, zero, w_max);
    W_hist_1pos_0charge[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_0q_%d", sec), Form("W_hist_1pos_0q_%d", sec), bins, zero, w_max);
    W_hist_1pos_gpart2[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_part2_%d", sec), Form("W_hist_1pos_part2_%d", sec), bins, zero, w_max);

    W_vs_q2_all_events[sec] =
        std::make_shared<TH2D>(Form("WQ2_sec_%d", sec), Form("WQ2_sec_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_%d", sec), Form("WQ2_1pos_%d", sec), bins, zero, w_max,
                                               bins, zero, q2_max);
    W_vs_q2_1pos_0charge[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_0q_%d", sec), Form("WQ2_1pos_0q_%d", sec), bins,
                                                       zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos_gpart2[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_part2_%d", sec), Form("WQ2_1pos_part2_%d", sec),
                                                      bins, zero, w_max, bins, zero, q2_max);

    W_hist_1pos_at90[sec] =
        std::make_shared<TH1D>(Form("W_1pos_at90_%d", sec), Form("W_1pos_at90_%d", sec), bins, zero, w_max);
    W_vs_q2_1pos_at90[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at90_%d", sec), Form("WQ2_1pos_at90_%d", sec), bins,
                                                    zero, w_max, bins, zero, q2_max);
  }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e) {
  MomVsBeta[0]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());
}

void Histogram::Write_MomVsBeta() {
    for (short p = 0; p < num_sectors; p++) {
    MomVsBeta[p]->SetXTitle("Momentum (GeV)");
    MomVsBeta[p]->SetYTitle("#beta");
    MomVsBeta[p]->SetOption("COLZ1");
    MomVsBeta[p]->Write();
  }
}
