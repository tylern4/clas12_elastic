/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");

  makeHists();
  Nsparce = std::make_shared<THnSparseD>("nsparce", "nsparce", 3, sparce_bins, sparce_xmin, sparce_xmax);
  makeHists_electron_cuts();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  Write_SF();
  Nsparce->Write();
  std::cout << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  Write_WvsQ2();

  std::cout << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();
  deltaT_proton[before_cut]->Write();
  deltaT_proton[after_cut]->Write();

  std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
  TDirectory* Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
  Electron_Cuts->cd();
  Write_Electron_cuts();

  std::cout << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::makeHists() {
  deltaT_proton[before_cut] = std::make_shared<TH2D>("DeltaTProton", "DeltaTProton", bins, zero, 10, bins, -5, 5);
  deltaT_proton[after_cut] =
      std::make_shared<TH2D>("DeltaTProton_cut", "DeltaTProton_cut", bins, zero, 10, bins, -5, 5);

  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    MissingMass[sec] =
        std::make_shared<TH1D>(Form("MM2_hist_sec_%d", sec), Form("MM2_hist_sec_%d", sec), bins, -w_max, w_max);

    mass_pi0_hist[before_cut][sec] =
        std::make_shared<TH1D>(Form("mass_pi0_hist_%d", sec), Form("mass_pi0_hist_%d", sec), bins, 0, 1.0);
    mass_pi0_hist[after_cut][sec] = std::make_shared<TH1D>(Form("mass_pi0_hist_aferPcuts_%d", sec),
                                                           Form("mass_pi0_hist_aferPcuts_%d", sec), bins, 0, 1.0);

    mass_pairCalc_hist[before_cut][sec] = std::make_shared<TH1D>(
        Form("mass_pairCalc_hist_%d", sec), Form("mass_pairCalc_hist_%d", sec), bins * 2, 0.0, 1.0);
    mass_pairCalc_hist[after_cut][sec] = std::make_shared<TH1D>(
        Form("mass_pairCalc_hist_aferPcuts_%d", sec), Form("mass_pairCalc_hist_aferPcuts_%d", sec), bins * 2, 0.0, 1.0);

    W_hist_all_events[sec] =
        std::make_shared<TH1D>(Form("W_hist_sec_%d", sec), Form("W_hist_sec_%d", sec), bins, zero, w_max);
    W_hist_1pos[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_%d", sec), Form("W_hist_1pos_%d", sec), bins, zero, w_max);
    W_hist_1pos_0charge[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_MM0_%d", sec), Form("W_hist_1pos_MM0_%d", sec), bins, zero, w_max);
    W_hist_1pos_noOther[sec] = std::make_shared<TH1D>(Form("W_hist_1pos_noOther_%d", sec),
                                                      Form("W_hist_1pos_noOther_%d", sec), bins, zero, w_max);

    W_vs_q2_all_events[sec] =
        std::make_shared<TH2D>(Form("WQ2_sec_%d", sec), Form("WQ2_sec_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_%d", sec), Form("WQ2_1pos_%d", sec), bins, zero, w_max,
                                               bins, zero, q2_max);
    W_vs_q2_1pos_0charge[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_MM0_%d", sec), Form("WQ2_1pos_MM0_%d", sec), bins,
                                                       zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos_noOther[sec] = std::make_shared<TH2D>(
        Form("WQ2_1pos_noOther_%d", sec), Form("WQ2_1pos_noOther_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    for (auto&& det : detector_name) {
      int d = detector_fill[det.first];

      ThetaVsP[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_pos_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_pos_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);

      ThetaVsP_lowW[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_lowW_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_lowW_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);

      ThetaVsPCalc[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_Calc_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_Calc_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);
      MomVsBeta[d][sec] =
          std::make_shared<TH2D>(Form("MomVsBeta_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsBeta_%s_%d", det.second.c_str(), sec), 500, zero, p_max, 500, zero, 1.2);
      Phie_vs_Phip[d][sec] =
          std::make_shared<TH2D>(Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec), 500, -PI, PI, 500, -PI, PI);
      Phie_Phip_hist[d][sec] =
          std::make_shared<TH1D>(Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec), 500, zero, 2 * PI);
      W_hist_1pos_at180[d][sec] =
          std::make_shared<TH1D>(Form("W_1pos_at180_%s_%d", det.second.c_str(), sec),
                                 Form("W_1pos_at180_%s_%d", det.second.c_str(), sec), bins, zero, w_max);
      W_vs_q2_1pos_at180[d][sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at180_%s_%d", det.second.c_str(), sec),
                                                          Form("WQ2_1pos_at180_%s_%d", det.second.c_str(), sec), bins,
                                                          zero, w_max, bins, zero, q2_max);
      W_hist_1pos_at180_MM[d][sec] =
          std::make_shared<TH1D>(Form("W_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                 Form("W_1pos_at180_MM_%s_%d", det.second.c_str(), sec), bins, zero, w_max);
      W_vs_q2_1pos_at180_MM[d][sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                                             Form("WQ2_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                                             bins, zero, w_max, bins, zero, q2_max);
    }
  }
}
void Histogram::makeHists_electron_cuts() {
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    auto type = cut.second.c_str();
    EC_sampling_fraction[c] =
        std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type), Form("EC_sampling_fraction%s", type), bins, p_min,
                               p_max, bins, zero, 1.0);
    vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type), Form("vz_position%s", type), bins, -40, 40);
    pcal_sec[c] =
        std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s", type), bins, -420, 420, bins, -420, 420);
    dcr1_sec[c] =
        std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s", type), bins, -180, 180, bins, -180, 180);
    dcr2_sec[c] =
        std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s", type), bins, -270, 270, bins, -270, 270);
    dcr3_sec[c] =
        std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s", type), bins, -320, 320, bins, -320, 320);
  }
}
void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d) {
  sf_hist->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}
void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[before_cut]->Fill(_d->vz(0));
  pcal_sec[before_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[before_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[before_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[before_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[before_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[after_cut]->Fill(_d->vz(0));
  pcal_sec[after_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[after_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[after_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[after_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[after_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}
void Histogram::Write_SF() { sf_hist->Write(); }
void Histogram::Write_Electron_cuts() {
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    vz_position[c]->SetXTitle("vz (GeV)");
    if (vz_position[c]->GetEntries()) vz_position[c]->Write();
    pcal_sec[c]->SetXTitle("x (cm)");
    pcal_sec[c]->SetYTitle("y (cm)");
    pcal_sec[c]->SetOption("COLZ1");
    if (pcal_sec[c]->GetEntries()) pcal_sec[c]->Write();

    dcr1_sec[c]->SetXTitle("x (cm)");
    dcr1_sec[c]->SetYTitle("y (cm)");
    dcr1_sec[c]->SetOption("COLZ1");
    if (dcr1_sec[c]->GetEntries()) dcr1_sec[c]->Write();

    dcr2_sec[c]->SetXTitle("x (cm)");
    dcr2_sec[c]->SetYTitle("y (cm)");
    dcr2_sec[c]->SetOption("COLZ1");
    if (dcr2_sec[c]->GetEntries()) dcr2_sec[c]->Write();

    dcr3_sec[c]->SetXTitle("x (cm)");
    dcr3_sec[c]->SetYTitle("y (cm)");
    dcr3_sec[c]->SetOption("COLZ1");
    if (dcr3_sec[c]->GetEntries()) dcr3_sec[c]->Write();

    EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
    EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
    EC_sampling_fraction[c]->SetOption("COLZ1");
    EC_sampling_fraction[c]->Write();
  }
}
void Histogram::Fill_Sparce(const std::shared_ptr<Reaction>& _e) {
  std::lock_guard<std::mutex> lk(mutex);
  double ret[NUM_DIM] = {_e->W(), _e->Q2(), static_cast<double>(_e->sec())};
  Nsparce->Fill(ret);
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  short pos_det = _e->pos_det();
  if ((sec > 0 && sec < NUM_SECTORS) || pos_det != -1) {
    W_hist_all_events[all_sectors]->Fill(_e->W());
    W_vs_q2_all_events[all_sectors]->Fill(_e->W(), _e->Q2());
    W_hist_all_events[sec]->Fill(_e->W());
    W_vs_q2_all_events[sec]->Fill(_e->W(), _e->Q2());

    if (_e->onePositive()) {
      W_hist_1pos[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos[sec]->Fill(_e->W());
      W_vs_q2_1pos[sec]->Fill(_e->W(), _e->Q2());

      Phie_vs_Phip[both_detectors][all_sectors]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[both_detectors][all_sectors]->Fill(_e->phi_diff());
      Phie_vs_Phip[both_detectors][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[both_detectors][sec]->Fill(_e->phi_diff());

      Phie_vs_Phip[pos_det][all_sectors]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][all_sectors]->Fill(_e->phi_diff());
      Phie_vs_Phip[pos_det][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][sec]->Fill(_e->phi_diff());
    }
    if (_e->onePositive_MM0()) {
      W_hist_1pos_0charge[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_0charge[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_0charge[sec]->Fill(_e->W());
      W_vs_q2_1pos_0charge[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_noOther()) {
      W_hist_1pos_noOther[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_noOther[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_noOther[sec]->Fill(_e->W());
      W_vs_q2_1pos_noOther[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180()) {
      MissingMass[all_sectors]->Fill(_e->MM2());
      MissingMass[sec]->Fill(_e->MM2());

      W_hist_1pos_at180[both_detectors][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180[both_detectors][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[both_detectors][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[both_detectors][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180[pos_det][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180_MM0()) {
      W_hist_1pos_at180_MM[both_detectors][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[both_detectors][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[both_detectors][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[both_detectors][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180_MM[pos_det][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
  }
}

void Histogram::Write_WvsQ2() {
  TDirectory* phi_folder = RootOutputFile->mkdir("Phi");
  phi_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      Phie_vs_Phip[j][i]->SetXTitle("Phie");
      Phie_vs_Phip[j][i]->SetYTitle("Phip");
      Phie_vs_Phip[j][i]->SetOption("COLZ");
      Phie_vs_Phip[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      Phie_Phip_hist[j][i]->SetXTitle("Phi");
      Phie_Phip_hist[j][i]->Write();
    }
  }
  phi_folder->Write();
  delete phi_folder;

  TDirectory* at180_folder = RootOutputFile->mkdir("at180");
  at180_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_hist_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_vs_q2_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180[j][i]->Write();
    }
  }
  at180_folder->Write();
  delete at180_folder;

  TDirectory* at180_MM_folder = RootOutputFile->mkdir("at180_MM");
  at180_MM_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_hist_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180_MM[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_vs_q2_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180_MM[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180_MM[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180_MM[j][i]->Write();
    }
  }
  at180_MM_folder->Write();
  delete at180_MM_folder;

  TDirectory* W_vs_Q2_folder = RootOutputFile->mkdir("W_vs_Q2");
  W_vs_Q2_folder->cd();
  for (int i = 0; i < NUM_SECTORS; i++) {
    MissingMass[i]->SetXTitle("MM^2 (GeV)");
    MissingMass[i]->Write();

    mass_pi0_hist[before_cut][i]->SetXTitle("Mass (GeV)");
    mass_pi0_hist[before_cut][i]->Write();

    mass_pi0_hist[after_cut][i]->SetXTitle("Mass (GeV)");
    mass_pi0_hist[after_cut][i]->Write();

    mass_pairCalc_hist[before_cut][i]->SetXTitle("Mass (GeV)");
    mass_pairCalc_hist[before_cut][i]->Write();

    mass_pairCalc_hist[after_cut][i]->SetXTitle("Mass (GeV)");
    mass_pairCalc_hist[after_cut][i]->Write();

    W_hist_all_events[i]->SetXTitle("W (GeV)");
    W_hist_all_events[i]->Write();
    W_hist_1pos[i]->SetXTitle("W (GeV)");
    W_hist_1pos[i]->Write();
    W_hist_1pos_0charge[i]->SetXTitle("W (GeV)");
    W_hist_1pos_0charge[i]->Write();
    W_hist_1pos_noOther[i]->SetXTitle("W (GeV)");
    W_hist_1pos_noOther[i]->Write();

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

    W_vs_q2_1pos_noOther[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_noOther[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_noOther[i]->SetOption("COLZ");
    W_vs_q2_1pos_noOther[i]->Write();
  }
  W_vs_Q2_folder->Write();
  delete W_vs_Q2_folder;
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e) {
  if (!_e->onePositive_at180()) return;
  if (_e->pos_det() == -1) return;

  MomVsBeta[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  MomVsBeta[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  if (_e->W() < 2.0) {
    ThetaVsP_lowW[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
  }
}

void Histogram::Write_MomVsBeta() {
  for (short i = 0; i < detector_name.size(); i++) {
    for (short p = 0; p < NUM_SECTORS; p++) {
      MomVsBeta[i][p]->SetXTitle("Momentum (GeV)");
      MomVsBeta[i][p]->SetYTitle("#beta");
      MomVsBeta[i][p]->SetOption("COLZ1");
      MomVsBeta[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsP[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP[i][p]->SetYTitle("#theta");
      ThetaVsP[i][p]->SetOption("COLZ1");
      ThetaVsP[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsP_lowW[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP_lowW[i][p]->SetYTitle("#theta");
      ThetaVsP_lowW[i][p]->SetOption("COLZ1");
      ThetaVsP_lowW[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsPCalc[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsPCalc[i][p]->SetYTitle("#theta");
      ThetaVsPCalc[i][p]->SetOption("COLZ1");
      ThetaVsPCalc[i][p]->Write();
    }
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt) {
  for (int i = 0; i < dt->gpart(); i++) {
    if (dt->charge(i) == 1) deltaT_proton[before_cut]->Fill(dt->mom(i), dt->dt_P(i));
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt, int part) {
  deltaT_proton[after_cut]->Fill(dt->mom(part), dt->dt_P(part));
}

void Histogram::Fill_pi0(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  short pos_det = _e->pos_det();

  for (auto&& m : _e->pair_mass()) {
    mass_pairCalc_hist[before_cut][all_sectors]->Fill(m);
    if ((sec > 0 && sec < NUM_SECTORS) || pos_det != -1) {
      mass_pairCalc_hist[before_cut][sec]->Fill(m);
      if (_e->onePositive_at180_MM0()) {
        mass_pairCalc_hist[after_cut][all_sectors]->Fill(m);
        mass_pairCalc_hist[after_cut][sec]->Fill(m);
      }
    }
  }

  if (_e->pi0_mass() < 0.0001) return;

  mass_pi0_hist[before_cut][all_sectors]->Fill(_e->pi0_mass());
  if ((sec > 0 && sec < NUM_SECTORS) || pos_det != -1) {
    mass_pi0_hist[before_cut][sec]->Fill(_e->pi0_mass());
    if (_e->onePositive_at180_MM0()) {
      mass_pi0_hist[after_cut][all_sectors]->Fill(_e->pi0_mass());
      mass_pi0_hist[after_cut][sec]->Fill(_e->pi0_mass());
    }
  }
}
