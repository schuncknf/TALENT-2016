#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "Solver.h"

/// class HartreeFock - 
class HartreeFock : public Solver {
  // Attributes
public:
  /// Matrix of eigenstates of the hamiltonian
  arma::field<arma::mat> D;
  /// Occupation numbers
  arma::field<arma::vec> occ;
  // Operations
public:
  HartreeFock (System & system);
  void calc (arma::field<arma::mat> & H);
  ~HartreeFock ();
};

#endif
