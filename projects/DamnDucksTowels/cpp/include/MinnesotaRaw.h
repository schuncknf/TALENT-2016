#ifndef MINNESOTARAW_H
#define MINNESOTARAW_H

#include <string>

#include "Interaction.h"
#include "FullSpBasis.h"

/// class MinnesotaRaw - Class describing the Minnesota pontential when considering all quantum numbers
class MinnesotaRaw : public Interaction
{
  // Attributes
public:
  /// Field of matrices containing the two-body matrix elements
  arma::field<arma::mat> TBME;
  /// Number of particle types considered
  int nParticleTypes;
  // Operations
public:
  /// Constructor for the MinnesotaRaw class
  MinnesotaRaw (FullSpBasis &_basis, int _nParticleTypes, std::string _filename);
  /// Destructor
  ~MinnesotaRaw ();
  /// Accesses the matrix elements of the interaction
  double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId);
  /// Returns a minimal set of information
  std::string info ();
  /// Returns a complete set of information
  std::string toString ();
};

#endif
