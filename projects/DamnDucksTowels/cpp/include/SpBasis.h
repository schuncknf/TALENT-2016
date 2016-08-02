#ifndef SPBASIS_H
#define SPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class SpBasis - Implements a harmonic oscillator single-particle basis
class SpBasis : public Basis
{
  // Attributes
public:
  /// Parameter of the harmonic oscillator basis
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed values for the quantum number l
  arma::ivec lMax;
  /// Maximum allowed values for the quantum number m
  arma::imat mMax;
  /// Auxiliary parameter
  double nu;
private:
  /// Normalization coefficients 
  arma::vec N;
  // Operations
public:
  /// Constructor of the SpBasis class
  SpBasis (double _omega, int _nMax, int _lMax);
  /// Destructor
  ~SpBasis ();
  /// Calculates the radial wave function
  void evalRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  /// Calculates the first derivative of the radial wave function
  void evalDerivativeRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  /// Alternate method for calculating the radial wave function
  void evalRadialWaveFunctionNoExp (arma::mat &wfMatrix, arma::vec &r);
  /// Kronecker delta for the spin
  int deltaSpin (int idx1, int idx2);
private:
  /// Calculating normalization coefficient
  void calcN ();
};

#endif
