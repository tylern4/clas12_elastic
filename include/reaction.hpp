/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::vector<std::unique_ptr<TLorentzVector>> _other;

  bool _hasE = false;
  bool _hasPos = false;
  bool _hasOther = false;

  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  int _total_charge = 0;
  short _pos_det = 0;
  short _sector = -1;

  float _pos_beta = NAN;

  float _MM = NAN;
  float _MM2 = NAN;

  float _W = NAN;
  float _Q2 = NAN;

  void SetElec();

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();

  void SetPositive(int i);
  void SetOther(int i);

  void CalcMissMass();
  float MM();
  float MM2();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }
  inline short sec() { return _data->dc_sec(0); }
  inline short det() { return abs(_data->status(0) / 1000); }
  inline float pos_beta() { return _pos_beta; }
  inline float pos_P() { return _prot->P(); }
  inline short pos_det() {
    if (_pos_det == 2) return 0;
    if (_pos_det == 4) return 1;
    return -1;
  }

  inline float phi_e() { return _elec->Phi(); }
  inline float phi_p() { return _prot->Phi(); }
  inline float phi_diff() { return abs(_elec->Phi() - _prot->Phi()); }
  inline bool phi_diff_90() {
    if (phi_diff() > 3.05 && phi_diff() < 3.2) return true;
    return false;
  }

  inline bool onePositive() { return (_hasE && _hasPos); }
  inline bool onePositive_MM0() { return (_hasE && _hasPos && abs(MM2()) < 0.1); }
  inline bool onePositive_part2() { return (_hasE && _hasPos && _data->gpart() == 2); }
  inline bool onePositive_at90() { return (_hasE && _hasPos && phi_diff_90()); }
  inline bool onePositive_at90_MM0() { return (_hasE && _hasPos && phi_diff_90() && abs(MM2()) < 0.1); }
};

#endif
