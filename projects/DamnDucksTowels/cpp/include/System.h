#include "Basis.h"

/// class System - 
class System {
  // Attributes
public:
  Basis* basis;
  /// Number of different particles
  arma::ivec particleNumbers ;
  /// Names of the different particle types
  std::vector<std::string> particleNames ;
  /// Type of potential used for calculations
  int potType ;
  /// Name of the system
  std::string name = "unnamed" ;

  // Methods
  /// First guess of the Hamiltonian
  virtual void calcH0 (arma::field<arma::mat> & H, int type) = 0;
  /// Calculation of the Hamiltonian
  virtual void calcH (arma::field<arma::mat> & H, arma::field<arma::mat> & RG) = 0;
};

