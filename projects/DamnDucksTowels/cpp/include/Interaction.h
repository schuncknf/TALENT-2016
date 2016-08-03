#ifndef INTERACTION_H
#define INTERACTION_H

#include <armadillo>

#include "Basis.h"

/// class Interaction - Abstract class describing the interaction
class Interaction
{
  // Attributes
public:
  /// Basis in which the interaction is computed
  arma::field<Basis *> particleBases;
  // Operations
public:
  /// Constructor for the Interaction class using multiple bases
  Interaction (arma::field<Basis *> _particleBases);
  /// Constructor for the Interaction class using a single basis with different particle types
  Interaction (Basis &_basis, int _nParticleTypes);
  /// Destructor
  virtual ~Interaction () = 0;
  /// Method accessing the matrix elements of the interaction
  double get (arma::field<arma::mat> &R, int ti, int i, int tj, int j);
  /// Alternate method accessing the matrix elements of the interaction
  double get (arma::field<arma::mat> &R, int ti0, int ti1, int i0, int i1, int tj0, int tj1, int j0, int j1);
  /// Abstract method accessing the matrix elements of the interaction
  virtual double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kid) = 0;
  /// Returns a minimal set of information about the object
  virtual std::string info () = 0;
  /// Returns a complete set of information about the object
  virtual std::string toString () = 0;
};

#endif
