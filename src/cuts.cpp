/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "cuts.hpp"
#include <iostream>

Cuts::Cuts(const std::shared_ptr<Branches12>& data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

bool Cuts::ElectronCuts() {
  bool _elec = true;
  // Number of good particles is greater than 0
  // So that we can check at(0) without errors
  _elec &= (_data->gpart() > 0);
  if (!_elec) return false;

  _elec &= (_data->charge(0) == NEGATIVE);
  _elec &= (_data->pid(0) == ELECTRON);
  _elec &= (_data->beta(0) > 0.05);
  //_elec &= !std::isnan(_data->cc_nphe_tot(0));
  _elec &= (_data->ec_tot_energy(0) / _data->p(0) > 0.2);
  _elec &= (_data->ec_tot_energy(0) / _data->p(0) < 0.3);

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
