#ifndef RAWINTERACTION_H
#define RAWINTERACTION_H

#include <armadillo>

#include "Interaction.h"
#include "Basis.h"

/// class RawInteraction - 
class RawInteraction : public Interaction {
  // Attributes
public:
  arma::field<arma::field<arma::field<arma::mat> > > potential;
  int nParticleTypes;
  // Operations
public:
  RawInteraction (Basis & _basis, int _nParticleTypes);
  ~RawInteraction ();
  double & set (int ti0, int i0, int ti1, int i1, int tj0, int j0, int tj1, int j1);
  double get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId);
};

#endif
