#ifndef MINNESOTAS0_H
#define MINNESOTAS0_H

#include <armadillo>
#include "ReducedSpBasis.h"

#include "Interaction.h"

/// class MinnesotaS0 - 
class MinnesotaS0 : public Interaction {
  // Attributes
public:
  arma::field<arma::mat> TBME;
  // Operations
public:
  MinnesotaS0 (ReducedSpBasis &_basis);
  ~MinnesotaS0 ();
  double get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId);
private:
  ReducedSpBasis &basis;
  void calculateTBME(int quadratureOrder);
  void printTBME();
  double V0R;
  double V0t;
  double V0s;
  double kR;
  double kt;
  double ks;
};

#endif
