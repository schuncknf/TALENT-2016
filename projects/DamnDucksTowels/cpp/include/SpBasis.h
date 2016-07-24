#ifndef SPBASIS_H
#define SPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class SpBasis - 
class SpBasis : public Basis {
  // Attributes
public:
  /// The harmonic oscillator frequency
  SpBasis(double _omega, int _nMax, int _lMax);
  ~SpBasis ();
  /// Parameter of the harmonic oscillator basis
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed value for the quantum number l
  arma::ivec lMax;
  /// Maximum allowed value for the quantum number m
  arma::imat mMax;
  // Operations
public:
  /// Calculating radial wave function
  void evalRadialWaveFunction (arma::mat &wfMatrix, arma::vec r);
  int deltaSpin (int idx1, int idx2);
private:
  /// Auxiliary parameter
  double nu;
  /// Calculating normalization coefficient
  void calcN();
  arma::vec N;
};

#endif
