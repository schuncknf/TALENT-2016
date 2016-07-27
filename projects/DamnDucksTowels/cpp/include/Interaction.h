#ifndef INTERACTION_H
#define INTERACTION_H

#include <armadillo>

#include "Basis.h"

/// class Interaction - 
class Interaction {
  // Attributes
public:
  arma::field<Basis *> particleBases;
  // Operations
public:
  Interaction (arma::field<Basis *> _particleBases);
  Interaction (Basis & _basis, int _nParticleTypes);
  virtual ~Interaction () = 0;
  double get (arma::field<arma::mat> & R, int ti, int i, int tj, int j);
  double get (arma::field<arma::mat> & R, int ti0, int ti1, int i0, int i1, int tj0, int tj1, int j0, int j1);
  virtual double get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kid) = 0;
  virtual std::string info () = 0;
  virtual std::string toString () = 0;
};

#endif
