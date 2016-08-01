#ifndef FULLSPBASIS_H
#define FULLSPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class FullSpBasis - 
class FullSpBasis : public Basis {
  // Attributes
public:
  /// Parameter of the harmonic oscillator basis
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed value for the quantum number l
  int lMax;
  /// Auxiliary parameter
  double nu;
private:
  arma::vec N;
  // Operations
public:
  FullSpBasis (double _omega, int _nMax, int _lMax);
  ~FullSpBasis ();
  /// Calculating radial wave function
  void evalRadialWaveFunction (arma::mat & wfMatrix, arma::vec & r);
  void evalDerivativeRadialWaveFunction (arma::mat & wfMatrix, arma::vec & r);
  void evalRadialWaveFunctionNoExp (arma::mat & wfMatrix, arma::vec & r);
private:
  /// Calculating normalization coefficient
  void calcN ();
};

#endif
