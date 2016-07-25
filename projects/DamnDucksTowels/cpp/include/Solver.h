#ifndef SOLVER_H
#define SOLVER_H

#include "System.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  System * system;
  /// Individual energies (indivEnergies(ptype, idState))
  arma::mat indivEnergies;
  double cvg;
  // Operations
public:
  Solver (System & _system, unsigned int _dNumber);
  virtual ~Solver () = 0;
  virtual void calc () = 0;
  std::string info ();
  std::string toString ();
};

#endif
