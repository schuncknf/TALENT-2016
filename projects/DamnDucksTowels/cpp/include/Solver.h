#ifndef SOLVER_H
#define SOLVER_H

#include "System.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  System * system;
  /// Deneralized densities: RG(dtype, ptype) where dtype is an integer representing 
  arma::field<arma::mat> RG;
  arma::vec indivEnergies;
  double cvg;
  // Operations
public:
  Solver (int matSize);
  virtual void calc ( arma::field<arma::mat> &) = 0;
  virtual ~Solver () = 0;
};

#endif
