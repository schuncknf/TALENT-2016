#include "System.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  System & system;
  arma::field<arma::mat> RG(2,system.particleNumbers.n_elem());
	// Operations
public:
  virtual void iter(arma::field<arma::mat>) = 0;
};

