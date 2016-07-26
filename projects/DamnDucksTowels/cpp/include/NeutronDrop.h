#ifndef NEUTRONDROP_H
#define NEUTRONDROP_H

#include <armadillo>

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class NeutronDrop - 
class NeutronDrop : public System {
  // Operations
public:
  NeutronDrop (Basis & _basis, Interaction & _inter, int _nbNeut);
  /// Destructor of the NeutronDrop class
  ~NeutronDrop ();
};

#endif
