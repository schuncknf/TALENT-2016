#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "System.h"

/// class Nucleus - 
class Nucleus : public System {
  // Operations
public:
  Nucleus (Basis & basis, int nbNeut, int nbProt);
  void calcH0 (arma::field<arma::mat> & H0, int type);
  void calcH (arma::field<arma::mat> & H, arma::field<arma::mat> & RG);
  void calcKineticField (arma::field<arma::mat> & RG);
  ~Nucleus ();
};

#endif
