/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;

class Histogram {
 protected:
  std::shared_ptr<TFile> RootOutputFile;
  std::shared_ptr<TCanvas> def;

  int bins = 500;
  double p_min = 0.0;
  double p_max = 10.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 8.0;
  double w_max = 4.0;

  double zero = 0.0;
  static const short num_sectors = 7;

  // Kinematics
  TH1D_ptr W_hist_all_events[num_sectors];
  TH1D_ptr W_hist_1pos[num_sectors];
  TH1D_ptr W_hist_1pos_0charge[num_sectors];
  TH1D_ptr W_hist_1pos_gpart2[num_sectors];
  TH1D_ptr W_hist_1pos_at90[num_sectors];

  TH2D_ptr W_vs_q2_all_events[num_sectors];
  TH2D_ptr W_vs_q2_1pos[num_sectors];
  TH2D_ptr W_vs_q2_1pos_0charge[num_sectors];
  TH2D_ptr W_vs_q2_1pos_gpart2[num_sectors];
  TH2D_ptr W_vs_q2_1pos_at90[num_sectors];

  TH2D_ptr MomVsBeta[num_sectors];

  TH2D_ptr Phie_vs_Phip[num_sectors];
  TH1D_ptr Phie_Phip_hist[num_sectors];

 public:
  Histogram(const std::string& output_file);
  ~Histogram();

  // W and Q^2
  void makeHists();
  void Fill_WvsQ2(const std::shared_ptr<Reaction>& _e);
  void Write_WvsQ2();

  // P and E
  void Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e);
  void Write_MomVsBeta();

  //
  void Write();
};

#endif
