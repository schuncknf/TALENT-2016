#ifndef NEUTRONDROP_H
#define NEUTRONDROP_H

#include <armadillo>

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class NeutronDrop - 
class NeutronDrop : public System {
  // Nested enumeration for initialization of the pot
  enum { SQUARE, WOOD_SAXON };
  // Operations
public:
  NeutronDrop (Basis & _basis, Interaction & _inter, int _nbNeut);
  /// Destructor of the NeutronDrop class
  ~NeutronDrop ();
  void calcH0 (int _type);
  /// calculate the hamiltonian of the system from densities
  void calcH ();
};

#endif
