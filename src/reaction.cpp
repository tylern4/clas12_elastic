/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;  // atof(getenv("CLAS12_E"));

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();
  _prot = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _total_charge += _data->charge(0);
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetPositive(int i) {
  _total_charge += _data->charge(i);
  _pos_beta = _data->beta(i);
  _numPos++;
  _hasPos = true;
  _pos_det = abs(_data->status(i) / 1000);
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}

void Reaction::SetOther(int i) {
  _total_charge += _data->charge(i);
  _numOther++;
  _hasOther = true;
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  *mm += (*_gamma + *_target);
  if (_prot != nullptr) {
    for (auto& _o : _other) *mm -= *_o;
    *mm -= *_prot;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}