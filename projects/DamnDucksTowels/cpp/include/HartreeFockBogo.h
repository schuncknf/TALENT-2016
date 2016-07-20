#include "Solver.h"

#include <armadillo>

/// class HartreeFockBogo - 
class HartreeFockBogo : public Solver {
  // Attributes
public:
  /// Density matrix
  arma::mat rho;
  /// Pairing matrix
  arma::mat kappa;
  /// Transformation matrix from HO to HF
  arma::mat D;
  /// Single-particle energies 
  arma::vec e;

// Operations
public:
  HartreeFockBogo(System & system);
  void iter(arma::field<arma::mat> H);
};

