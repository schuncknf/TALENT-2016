#include "System.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  System * system;
  /// Emulation of generalized density matrix
  arma::field<arma::mat> RG;
  /// Vector for the single-particle energies
  arma::vec indivEnergies;
	// Operations
public:
  virtual void iter(arma::field<arma::mat>) = 0;
};

