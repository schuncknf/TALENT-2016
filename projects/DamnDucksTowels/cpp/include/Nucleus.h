#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "System.h"
#include "Interaction.h"
#include "Basis.h"

/// class Nucleus - 
class Nucleus : public System {
public:
  int nPoints;
  // Operations
public:
  Nucleus (Basis & _basis, Interaction & _inter, int _nbNeut, int _nbProt, int _nPoints = 50);
  ~Nucleus ();
  void calcH ();
};

#endif
