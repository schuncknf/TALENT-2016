#ifndef SOLVER_H
#define SOLVER_H

#include "System.h"

/// class Solver - Abstract class representing a generic solver
class Solver
{
  // Attributes
public:
  System *system;
  /// Individual energies (indivEnergies(ptype, idState))
  arma::mat indivEnergies;
  /// Maximum acceptable difference in individual energy between two iterations at convergence
  double cvg;
  // Operations
public:
  /// Constructor of the Solver class
  Solver (System &_system);
  /// Destructor
  virtual ~Solver () = 0;
  /// Run a sole iteration of the solver
  virtual void run () = 0;
  /// Initialize the Hamiltonian H
  void initH (int type);
  /// Compute the new value of the Hamiltonian H
  virtual void calcH () = 0;
  /// Returns a minimum set of information
  std::string info ();
  /// Returns a full set of information
  std::string toString ();
};

#endif
