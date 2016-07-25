#ifndef MINNESOTAS0_H
#define MINNESOTAS0_H

#include <armadillo>

#include "Interaction.h"

/// class MinnesotaS0 - 
class MinnesotaS0 : public Interaction {
  // Attributes
public:
  arma::field<arma::field<arma::field<arma::mat> > > TBME;
  // Operations
public:
  MinnesotaS0 ();
  ~MinnesotaS0 ();
  double get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId);
};

#endif
