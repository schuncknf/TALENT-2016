#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>
#include <armadillo>

#include "Basis.h"

/// class System - 
class System {
  // Attributes
public:
  /// The name of the system
  std::string name;
  Basis * basis;
  /// Number of particles for each particles type (i.e. neutrons, protons, nucleus, e
  arma::ivec particleNumbers;
  /// Name of each particle type (i.e. "neutron", "proton", "nucleus", "electron", et
  std::vector<std::string> particleNames;
  arma::field<arma::mat> TBME;
  // Operations
public:
  System (std::string _name, Basis & _basis, arma::ivec _particleNumbers, std::vector<std::string> particleNames, arma::field<arma::mat> & _TBME);
  virtual void calcH0 (arma::field<arma::mat> & H0, int type) = 0;
  virtual void calcH (arma::field<arma::mat> & H, arma::field<arma::mat> & RG) = 0;
  virtual ~System () = 0;
};

#endif
