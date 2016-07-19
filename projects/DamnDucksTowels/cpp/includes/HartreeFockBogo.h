#include "Solver.h"

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
  /// Occupation numbers
  arma::vec occ;
};

