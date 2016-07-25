#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class Nucleus - 
class Nucleus : public System {
  // Operations
public:
  Nucleus (Basis & _basis, Interaction & _inter, int _nbNeut, int _nbProt);
  ~Nucleus ();
  void calcH0 (int _type);
  void calcH ();
};

#endif
