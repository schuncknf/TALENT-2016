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
  /// Occupation numbers
  arma::vec occ;
};

