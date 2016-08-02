#ifndef NEUTRONDROP_H
#define NEUTRONDROP_H

#include <armadillo>

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class NeutronDrop - Implements a neutron drop system
class NeutronDrop : public System
{
  // Attributes
public:
  /// Omega parameter for the harmonic part of the Hamiltonian 
  double omega;
  // Operations
public:
  /// Constructor of the NeutronDrop class
  NeutronDrop (Basis &_basis, Interaction &_inter, int _nbNeut, double _omega, int _nPoints = 50);
  /// Destructor of the NeutronDrop class
  ~NeutronDrop ();
  /// Computes the Hamiltonian matrix
  void calcH ();
};

#endif
