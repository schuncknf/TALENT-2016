#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class Nucleus - Implements a nucleus system (Empty for the moment)
class Nucleus : public System
{
  // Operations
public:
  /// Constructor of the Nucleus class
  Nucleus (Basis &_basis, Interaction &_inter, int _nbNeut, int _nbProt, int _nPoints = 50);
  /// Destructor
  ~Nucleus ();
  /// Computes the Hamiltonian matrix
  void calcH ();
};

#endif
