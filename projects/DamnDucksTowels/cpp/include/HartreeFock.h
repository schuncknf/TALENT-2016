#include "Solver.h"
#include <armadillo>
/// class HartreeFock - 
class HartreeFock : public Solver {
  // Attributes
public:
  /// Matrix of eigenstates of the hamiltonian
  arma::mat D;
  /// Vector with the occupation numbers of the states
  arma::vec occ;

// Operations
public:
  HartreeFock(System & system);
  ~HartreeFock();
  void iter(arma::field<arma::mat> H);
};

