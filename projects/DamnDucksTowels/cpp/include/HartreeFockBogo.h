#ifndef HARTREEFOCKBOGO_h
#define HARTREEFOCKBOGO_h

#include <armadillo>

#include "Solver.h"
#include "System.h"
#include "FullSpBasis.h"
#include "MinnesotaRaw.h"

/// class HartreeFockBogo - Class representing a Hartree-Fock-Bogoliubov solver
class HartreeFockBogo : public Solver
{
  // Attributes
public:
  /// Transformation matrix from HO to HF
  arma::field<arma::mat> D;
  /// Occupation numbers
  arma::field<arma::vec> occ;
  // Operations
public:
  /// Constructor for the Hartree-Fock-Bogoliubov solver
  HartreeFockBogo (System &system);
  /// Destructor
  ~HartreeFockBogo ();
  /// Run a sole iteration of the Hartree-Fock-Bogoliubov solver
  void run ();
  /// Compute the Hamiltonian H (to be suppressed)
  void calcH ();
};

#endif
