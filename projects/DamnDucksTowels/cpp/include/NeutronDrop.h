#ifndef NEUTRONDROP_H
#define NEUTRONDROP_H

#include <armadillo>

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class NeutronDrop -
class NeutronDrop : public System
{
  // Attributes
public:
  double omega;
  // Operations
public:
  NeutronDrop (Basis &_basis, Interaction &_inter, int _nbNeut, double _omega, int _nPoints = 50);
  /// Destructor of the NeutronDrop class
  ~NeutronDrop ();
  void calcH ();
};

#endif
