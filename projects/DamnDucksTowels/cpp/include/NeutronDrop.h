#ifndef NEUTRONDROP_H
#define NEUTRONDROP_H

#include "System.h"
#include "Basis.h"

/// class NeutronDrop - 
class NeutronDrop : public System {
  // Nested enumeration for initialization of the pot
  enum { SQUARE, WOOD_SAXON };
  // Operations
public:
  NeutronDrop (int nbNeut, Basis & basis);
  void calcH0 (arma::field<arma::mat> & H0, int type);
  void calcH (arma::field<arma::mat> & H, arma::field<arma::mat> & RG);
  ~NeutronDrop ();
private:
  void calcKineticField (arma::field<arma::mat> & RG);
};

#endif
