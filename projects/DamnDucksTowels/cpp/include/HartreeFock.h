#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "Solver.h"
#include "System.h"
#include "quadrature.h"
#include "FullSpBasis.h"
#include "MinnesotaRaw.h"
#include "global.h"

/// class HartreeFock - Class representing a Hartree-Fock solver
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
  /// Constructor for the Hartree-Fock solver
  HartreeFock (System &system);
  /// Destructor
  ~HartreeFock ();
  /// Run a sole iteration of the Hartree-Fock solver
  void run ();
  /// Compute the Hamiltonian H (to be suppressed)
  void calcH ();
};

#endif
