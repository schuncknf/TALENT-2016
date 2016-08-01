#ifndef MINNESOTARAW_H
#define MINNESOTARAW_H

#include <string>

#include "Interaction.h"
#include "FullSpBasis.h"

/// class MinnesotaRaw -
class MinnesotaRaw : public Interaction
{
  // Attributes
public:
  arma::field<arma::mat> TBME;
  int nParticleTypes;
  // Operations
public:
  MinnesotaRaw (FullSpBasis &_basis, int _nParticleTypes, std::string _filename);
  ~MinnesotaRaw ();
  double get (arma::field<arma::mat> &R, arma::ivec &bType, arma::ivec &bId, arma::ivec &kType, arma::ivec &kId);
  std::string info ();
  std::string toString ();
};

#endif
