#include "Solver.h"

#include <armadillo>

/// class HartreeFockBogo - 
class HartreeFockBogo : public Solver {
  // Attributes
public:
  /// Transformation matrix from HO to HF
  arma::mat D;
  /// Vector with the occupation numebrs of the states 
  arma::vec occ;

// Operations
public:
  HartreeFockBogo(System & system);
  ~HartreeFockBogo();
  void iter(arma::field<arma::mat> H);
};

