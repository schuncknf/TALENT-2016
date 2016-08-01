#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "Solver.h"
#include "System.h"
#include "quadrature.h"
#include "FullSpBasis.h"
#include "MinnesotaRaw.h"
#include "global.h"

/// class HartreeFock -
class HartreeFock : public Solver
{
  // Attributes
public:
  /// Matrix of eigenstates of the hamiltonian
  arma::field<arma::mat> D;
  /// Occupation numbers
  arma::field<arma::vec> occ;
  // Operations
public:
  HartreeFock (System &system);
  ~HartreeFock ();
  void run ();
  void calcH ();
};

#endif
