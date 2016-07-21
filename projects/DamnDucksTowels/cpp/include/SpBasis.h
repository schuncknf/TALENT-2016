#ifndef SPBASIS_H
#define SPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class SpBasis - 
class SpBasis : public Basis {
  // Attributes
public:
  /// The harmonic oscillator frequency
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed value for the quantum number l
  arma::ivec lMax;
  arma::imat mMax;
  arma::vec radialWaveFunction;
private:
  double nu;
  arma::mat N;
  // Operations
public:
  SpBasis (double omega, int nMax, int lMax);
  void evalRadialWaveFunction (arma::mat & wfMatrix, arma::vec r);
  ~SpBasis ();
private:
  void calcN ();
};

#endif
