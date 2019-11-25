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
  std::unique_ptr<TLorentzVector> _prot = nullptr;
  std::vector<std::unique_ptr<TLorentzVector>> _pos;
  std::vector<std::unique_ptr<TLorentzVector>> _other;

  bool _hasE = false;
  bool _hasPos = false;
  bool _hasOther = false;

  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  std::vector<float> _pos_beta;
  std::vector<short> _pos_det;
  short _sector = -1;

  float _MM = NAN;
  float _MM2 = NAN;

  float _W = NAN;
  float _Q2 = NAN;

  void SetElec();

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12>& data, float beam_energy);
  ~Reaction();

  bool PosStats();
  void SetPositive(int i);
  void SetOther(int i);

  void CalcMissMass();
  float MM();
  float MM2();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }
  inline short sec() { return _data->dc_sec(0); }
  inline short det() { return abs(_data->status(0) / 1000); }
  inline float pos_beta() {
    if (_pos_beta.size() == 0) return NAN;
    return _pos_beta.front();
  }
  inline float pos_P() {
    if (_pos.size() == 0) return NAN;
    return _pos.front()->P();
  }
  inline short pos_det() {
    if (_pos_det.size() == 0) return -1;
    return detector_fill[_pos_det.front()];
  }
  inline float pos_theta() {
    if (_pos.size() == 0) return NAN;
    return _pos.front()->Theta();
  }

  inline float phi_e() { return _elec->Phi(); }
  inline float phi_p() {
    if (_pos.size() == 0) return NAN;
    return _pos.front()->Phi();
  }
  inline float phi_diff() { return abs(_elec->Phi() - phi_p()); }
  inline bool phi_diff_180() {
    // Cut around 10% of peak
    if (phi_diff() > (PI * 0.9) && phi_diff() < (PI * 1.1)) return true;
    return false;
  }
  inline bool MM_cut() { return abs(MM2()) < 0.1; }

  inline bool onePositive() { return (_hasE && _hasPos); }
  inline bool onePositive_MM0() { return (onePositive() && MM_cut()); }
  inline bool onePositive_part2() { return (onePositive() && _pos.size() == 1); }
  inline bool onePositive_at180() { return (onePositive() && phi_diff_180()); }
  inline bool onePositive_at180_MM0() { return (onePositive_at180() && MM_cut()); }
};

#endif
