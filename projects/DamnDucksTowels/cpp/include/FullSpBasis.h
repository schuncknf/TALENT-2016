#ifndef FULLSPBASIS_H
#define FULLSPBASIS_H

#include <armadillo>

#include "Basis.h"

/// class FullSpBasis - Describes a single-particle basis with all quantum numbers
class FullSpBasis : public Basis
{
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
  /// Normalization coefficients for the basis states 
  arma::vec N;
  // Operations
public:
  /// Constructor for the FullSpBasis class
  FullSpBasis (double _omega, int _nMax, int _lMax);
  /// Destructor
  ~FullSpBasis ();
  /// Calculates the radial wave function
  void evalRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  /// Calculates the first derivative of the radial wave expansion
  void evalDerivativeRadialWaveFunction (arma::mat &wfMatrix, arma::vec &r);
  /// Calculates the radial wave function with an alternate method
  void evalRadialWaveFunctionNoExp (arma::mat &wfMatrix, arma::vec &r);
private:
  /// Calculates normalization coefficients
  void calcN ();
};

#endif
