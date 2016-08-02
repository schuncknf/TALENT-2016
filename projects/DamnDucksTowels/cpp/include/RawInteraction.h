#ifndef RAWINTERACTION_H
#define RAWINTERACTION_H

#include <armadillo>

#include "Interaction.h"
#include "Basis.h"

/// class RawInteraction -
class RawInteraction : public Interaction
{
  // Attributes
public:
  /// Field of field of matrices to store the matrix elements of the pontential
  arma::field<arma::field<arma::field<arma::mat> > > potential;
  /// Number of different particle types considered
  int nParticleTypes;
  // Operations
public:
  /// Constructor of the RawInteraction class
  RawInteraction (Basis &_basis, int _nParticleTypes);
  /// Destructor
  ~RawInteraction ();
  /// Set the matrix elements of the interaction
  double &set (int ti0, int i0, int ti1, int i1, int tj0, int j0, int tj1, int j1);
  /// Accesses the matrix element of the interaction
  double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId);
  /// Returns a minimal set of information
  std::string info ();
  /// Returns a complete set of information
  std::string toString ();
};

#endif
