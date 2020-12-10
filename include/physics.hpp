/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include <TLorentzVector.h>
#include <map>
#include <unordered_map>
#include "TROOT.h"
#include "constants.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
double xb_calc(double Q2, double E_prime);
// overload with 4 vectors instaed of other calculations
double xb_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
float theta_calc_rad(const float &cosz);
float phi_calc_rad(const float &cosx, const float &cosy);

double vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta);
double deltat(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r);
}  // namespace physics

#endif
