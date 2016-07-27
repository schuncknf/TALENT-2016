#ifndef MINNESOTAS0_H
#define MINNESOTAS0_H

#include <armadillo>

#include "Interaction.h"
#include "ReducedSpBasis.h"

/// class MinnesotaS0 - 
class MinnesotaS0 : public Interaction {
  // Attributes
public:
  arma::field<arma::mat> TBME;
  int nPoints;
private:
  ReducedSpBasis & basis;
  // Operations
public:
  MinnesotaS0 (ReducedSpBasis & _basis, int _nPoints = 120);
  ~MinnesotaS0 ();
  double get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId);
  void calc ();
  std::string info ();
  std::string toString ();
};

#endif
