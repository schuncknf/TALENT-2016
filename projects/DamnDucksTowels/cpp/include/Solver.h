#ifndef SOLVER_H
#define SOLVER_H

#include "System.h"

/// class Solver -
class Solver
{
  // Attributes
public:
  System *system;
  /// Individual energies (indivEnergies(ptype, idState))
  arma::mat indivEnergies;
  double cvg;
  // Operations
public:
  Solver (System &_system);
  virtual ~Solver () = 0;
  virtual void run () = 0;
  void initH (int type);
  virtual void calcH () = 0;
  std::string info ();
  std::string toString ();
};

#endif
