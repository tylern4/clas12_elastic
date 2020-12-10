/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"
#include "Math/VectorUtil.h"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, const float& beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetPositive(int i) {
  _pos_beta.push_back(_data->beta(i));
  _numPos++;
  _hasPos = true;
  _pos_det.push_back(abs(_data->status(i) / 1000));
  _pos.push_back(std::make_unique<TLorentzVector>());
  _pos.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}

bool Reaction::PosStats() {
  if (_pos.size() != 2) return false;
  if (abs(_pos.front()->Phi() - _pos.back()->Phi()) > 0.5) return false;
  if (abs(_pos.front()->Theta() - _pos.back()->Theta()) > 0.5) return false;
  if (_pos_det.front() == _pos_det.back()) return false;
  /*
  std::cout << "num_pos: " << _pos.size() << std::endl;
  std::cout << "num_phi: ";
  for (auto& _d : _pos_det) std::cout << detector_name[_d] << "\t";
  std::cout << std::endl;
  std::cout << "pos_phi: ";
  for (auto& _p : _pos) std::cout << _p->Phi() << "\t";
  std::cout << std::endl;
  std::cout << "pos_theta: ";
  for (auto& _p : _pos) std::cout << _p->Theta() << "\t";
  std::cout << std::endl;
  */

  return true;
}

void Reaction::SetOther(int i) {
  _numOther++;
  _hasOther = true;
  if (_data->charge(i) != 0) {
    _other.push_back(std::make_unique<TLorentzVector>());
    _other.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
  if (_data->pid(i) == PHOTON) {
    _photons.push_back(std::make_unique<TLorentzVector>());
    _photons.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), 0);
    //_photons.back()->SetPxPyPzE(_data->px(i), _data->py(i), _data->pz(i), _data->ec_tot_energy(i));
  }
}

void Reaction::CalcMissMass() {
  if (_pos.size() > 0) {
    auto mm = std::make_unique<TLorentzVector>();
    *mm += (*_gamma + *_target);
    // for (auto& _p : _photons) *mm -= *_p;
    for (auto& _p : _pos) *mm -= *_p;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

float Reaction::MM() {
  if (isnan(_MM)) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (isnan(_MM2)) CalcMissMass();
  return _MM2;
}

void Reaction::CalcMassPi0() {
  _pi0_mass = 0;
  if (_photons.size() > 1) {
    auto mass = std::make_unique<TLorentzVector>();
    for (auto& _p : _photons) *mass += *_p;
    _pi0_mass = mass->M();
  }
}

float Reaction::pi0_mass() {
  if (isnan(_pi0_mass)) CalcMassPi0();
  return _pi0_mass;
}

void Reaction::CalcMassPairs() {
  std::lock_guard<std::mutex> lk(mutex);
  float _min_photon_E = 1.5;
  if (_photons.size() >= 2) {
    // Reverse photon vector
    std::vector<std::shared_ptr<TLorentzVector>> r_photons(_photons.rbegin(), _photons.rend());
    // For each photon
    for (auto&& _p : _photons) {
      // Remove last photon in reversed vector aka first photon now _p
      r_photons.pop_back();
      // For each photn in revered list
      for (auto& _rp : r_photons) {
        // Cut on minimum photon energy
        if (_p->Energy() < _min_photon_E || _rp->Energy() < _min_photon_E) continue;
        auto phi = _p->Angle(_rp->Vect());
        if (phi < 0.1) continue;
        // Add together pair and get mass
        _pair_mass.push_back((*_p + *_rp).M());
      }
    }
  }
}

std::vector<float> Reaction::pair_mass() {
  if (_pair_mass.size() == 0) CalcMassPairs();
  return _pair_mass;
}