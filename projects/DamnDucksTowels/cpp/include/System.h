#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>
#include <armadillo>

#include "Interaction.h"
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
  Interaction * inter;
  /// Hamiltonian of the system used as H(dtype, ptype) -> cf R comment
  arma::field<arma::mat> H;
  /// densities of the system, used as R(dtype, ptype) where dtype is the density typ
  arma::field<arma::mat> R;
  // Operations
public:
  System (std::string _name, Basis & _basis, arma::ivec _particleNumbers, std::vector<std::string> _particleNames, Interaction & _inter);
  virtual ~System () = 0;
  arma::mat & getH (int dType, int pType);
  arma::mat & getR (int dType, int pType);
  std::string info ();
  std::string toString ();
protected:
  void addParticleType (std::string name, int number);
};

#endif
