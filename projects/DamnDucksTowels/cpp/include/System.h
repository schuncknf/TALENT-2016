#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>
#include <armadillo>

#include "Interaction.h"
#include "Basis.h"

/// class System - Implements the considered system
class System
{
  // Attributes
public:
  /// Name of the system
  std::string name;
  /// Basis of interest
  Basis *basis;
  /// Number of particles for each particles type (i.e. neutrons, protons, nucleus, e
  arma::ivec particleNumbers;
  /// Name of each particle type (i.e. "neutron", "proton", "nucleus", "electron", et
  std::vector<std::string> particleNames;
  /// Interaction of interest
  Interaction *inter;
  // Hamiltonian of the system used as H(dtype, ptype) -> cf R comment
  /// Field of matrices for the kinetic part of the Hamiltonian
  arma::field<arma::mat> Kinetic;
  /// Field of matrices for the two-body part of the Hamiltonian
  arma::field<arma::mat> Gamma;
  /// Field of matrices for the pairing part of the Hamiltonian
  arma::field<arma::mat> Delta;
  /// Densities of the system, used as R(dtype, ptype) where dtype is the density type
  arma::field<arma::mat> R;
  /// Number of points used for the quadrature
  int nPoints;
  // Operations
public:
  /// Constructor of the system class
  System (std::string _name, Basis &_basis, arma::ivec _particleNumbers, std::vector<std::string> _particleNames, Interaction &_inter, int _nPoints = 50);
  /// Destructor
  virtual ~System () = 0;
  /// Virtual mehtod to compute the Hamiltonian matrix
  virtual void calcH () = 0;
  /// Accesses the Kinetic part of the Hamiltonian for the pType particle type
  arma::mat &getKinetic (int pType);
  /// Accesses the Gamma part of the Hamiltonian for the pType particle type
  arma::mat &getGamma (int pType);
  /// Accesses the Delta (pairing) part of the Hamiltonian for the pType particle type
  arma::mat &getDelta (int pType);
  /// Accesses the different density matrices of the system (rho, kappa) for the pType particle type
  arma::mat &getR (int dType, int pType);
  /// Returns a minimal set of information
  std::string info ();
  /// Returns a complete set of information
  std::string toString ();
protected:
  /// Add a new type of particle named inside "name" and identified through "number"
  void addParticleType (std::string name, int number);
};

#endif
