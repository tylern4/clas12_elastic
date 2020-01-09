/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include <iostream>
#include "cuts.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12>& data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

bool Cuts::ElectronCuts() {
  bool _elec = true;
  // Number of good particles is greater than 0
  // So that we can check at(0) without errors
  _elec &= (_data->gpart() > 0);
  if (!_elec) return false;

  _elec &= (_data->gpart() < 20);

  _elec &= (_data->charge(0) == NEGATIVE);
  _elec &= (_data->pid(0) == ELECTRON);
  //_elec &= (_data->beta(0) > 0.05);
  _elec &= (_data->p(0) > 1.0);
  _elec &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  _elec &= (_data->vz(0) > -7.9 && _data->vz(0) < 2.0);
  _elec &= !std::isnan(_data->cc_nphe_tot(0));
  _elec &= (_data->ec_tot_energy(0) / _data->p(0) > 0.2);
  _elec &= (_data->ec_tot_energy(0) / _data->p(0) < 0.3);
  _elec &= (abs(_data->chi2pid(0)) < 3);

  float x_PCAL_rot, y_PCAL_rot, angle, height_PCAL, slope_PCAL, left_PCAL, right_PCAL, radius2_PCAL, dcR1_height,
      dcR2_height, dcR3_height, x1_rot, y1_rot, x2_rot, y2_rot, x3_rot, y3_rot, slope, left_r1, right_r1, left_r2,
      right_r2, left_r3, right_r3, radius2_DCr1, radius2_DCr2, radius2_DCr3;

  x_PCAL_rot = _data->ec_pcal_y(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180) +
               _data->ec_pcal_x(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  y_PCAL_rot = _data->ec_pcal_y(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180) -
               _data->ec_pcal_x(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  angle = 60;
  height_PCAL = 45;  // stafen 45
  slope_PCAL = 1 / tan(0.5 * angle * PI / 180);
  left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
  right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
  radius2_PCAL = pow(height_PCAL + 6, 2) - pow(y_PCAL_rot, 2);  // circle radius r^2 = x^2 +
                                                                // y^2
  _elec &= (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && pow(x_PCAL_rot, 2) > radius2_PCAL && x_PCAL_rot < 372);
  dcR1_height = 25;  // 31 //20
  dcR2_height = 42;  // 47 // 30
  dcR3_height = 49;  // 53 // 40

  x1_rot = _data->dc_r1_y(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180) +
           _data->dc_r1_x(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  y1_rot = _data->dc_r1_y(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180) -
           _data->dc_r1_x(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  x2_rot = _data->dc_r2_y(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180) +
           _data->dc_r2_x(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  y2_rot = _data->dc_r2_y(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180) -
           _data->dc_r2_x(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  x3_rot = _data->dc_r3_y(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180) +
           _data->dc_r3_x(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180);
  y3_rot = _data->dc_r3_y(0) * cos((_data->dc_sec(0) - 1) * 60.0 * PI / 180) -
           _data->dc_r3_x(0) * sin((_data->dc_sec(0) - 1) * 60.0 * PI / 180);

  slope = 1 / tan(0.5 * angle * PI / 180);

  left_r1 = (dcR1_height - slope * y1_rot);
  right_r1 = (dcR1_height + slope * y1_rot);
  left_r2 = (dcR2_height - slope * y2_rot);
  right_r2 = (dcR2_height + slope * y2_rot);
  left_r3 = (dcR3_height - slope * y3_rot);
  right_r3 = (dcR3_height + slope * y3_rot);

  radius2_DCr1 = pow(32, 2) - pow(y1_rot, 2);  // 32 stafen // 21
  radius2_DCr2 = pow(49, 2) - pow(y2_rot, 2);  // 49 // 30
  radius2_DCr3 = pow(54, 2) - pow(y3_rot, 2);  // 54  // 40
  //
  _elec &= (x1_rot > left_r1 && x1_rot > right_r1 && pow(x1_rot, 2) > radius2_DCr1);
  _elec &= (x2_rot > left_r2 && x2_rot > right_r2 && pow(x2_rot, 2) > radius2_DCr2);
  _elec &= (x3_rot > left_r3 && x3_rot > right_r3 && pow(x3_rot, 2) > radius2_DCr3);
  return _elec;
}

bool Cuts::IsPip(int i) {
  if (_data->gpart() <= i) return false;
  bool _pip = true;
  _pip &= (_data->charge(i) == POSITIVE);
  _pip &= !(abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.2);
  return _pip;
}

bool Cuts::IsProton(int i) {
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  _proton &= (_data->charge(i) == POSITIVE);
  //_proton &= _data->pid(i) == PROTON;
  //_proton = (_data->p(i) > 1.0);
  if (!std::isnan(_dt->dt_P(i))) _proton &= (abs(_dt->dt_P(i)) < 0.5);
  // if (!std::isnan(_dt->dt_ctof_P(i))) _proton &= (abs(_dt->dt_ctof_P(i)) < 0.2);
  return _proton;
}

bool Cuts::IsPim(int i) {
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  return _pim;
}
