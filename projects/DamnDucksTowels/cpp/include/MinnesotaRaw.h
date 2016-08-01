#ifndef MINNESOTARAW_H
#define MINNESOTARAW_H

#include "Interaction.h"

/// class MinnesotaRaw -
class MinnesotaRaw : public MinnesotaS0
{
  // Attributes
public:
  arma::field<arma::field<arma::field<arma::mat> > > potential;
  int nParticleTypes;
  // Operations
public:
  MinnesotaRaw (FullSpBasis &_basis, int _nParticleTypes);
  ~MinnesotaRaw ();
  double &set (arma::vec &pTypes, arma::vec &n, arma::vec &l, arma::vec &j, arma::vec &mj);
  double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId);
  std::string info ();
  std::string toString ();
};

#endif
