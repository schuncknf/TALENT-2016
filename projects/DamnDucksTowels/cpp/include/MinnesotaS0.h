#ifndef MINNESOTAS0_H
#define MINNESOTAS0_H

#include <armadillo>

#include "Interaction.h"
#include "ReducedSpBasis.h"

/// class MinnesotaS0 - Class implementing the Minnesota potential when considering only the S-wave
class MinnesotaS0 : public Interaction
{
  // Attributes
public:
  /// Field of matrices containing the two-body matrix elements
  arma::field<arma::mat> TBME;
  /// Number of points used for the quadrature operation
  int nPoints;
private:
  /// Reduced single-particle basis considered
  ReducedSpBasis &basis;
  // Operations
public:
  /// Constructor for the MinnesotaS0 class
  MinnesotaS0 (ReducedSpBasis &_basis, int _nPoints = 120);
  /// Destructor
  ~MinnesotaS0 ();
  /// Accesses the matrix elements of the interaction
  double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId);
  /// Calculates the two-body matrix elements of the Minnesota potential considering only the S-wave
  void calc ();
  /// Returns a minimal set of information
  std::string info ();
  /// Returns a complete set of information
  std::string toString ();
};

#endif
