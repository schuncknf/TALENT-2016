#ifndef SPBASIS_H
#define SPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class SpBasis -
class SpBasis : public Basis
{
  // Attributes
public:
  /// Parameter of the harmonic oscillator basis
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed value for the quantum number l
  arma::ivec lMax;
  /// Maximum allowed value for the quantum number m
  arma::imat mMax;
  /// Auxiliary parameter
  double nu;
private:
  arma::vec N;
  // Operations
public:
  SpBasis (double _omega, int _nMax, int _lMax);
  ~SpBasis ();
  /// Calculating radial wave function
  void evalRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  void evalDerivativeRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  void evalRadialWaveFunctionNoExp (arma::mat &wfMatrix, arma::vec &r);
  int deltaSpin (int idx1, int idx2);
private:
  /// Calculating normalization coefficient
  void calcN ();
};

#endif
