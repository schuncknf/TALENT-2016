#include "Solver.h"
#include <armadillo>
/// class HartreeFock - 
class HartreeFock : public Solver {
  // Attributes
public:
  /// Density matrix
  arma::mat rho;
  /// Matrix of eigenstates of the hamiltonian
  arma::mat D;
  /// Single-particle energies 
  arma::vec e;

// Operations
public:
  HartreeFock(System & system);
  void iter(arma::field<arma::mat> H);
};

