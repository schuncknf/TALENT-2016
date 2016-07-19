#include "Basis.h"

/// class SpBasis - 
class SpBasis : public Basis {
  // Attributes
public:
  double omega;
  /// Maximum allowed value for the quantum number n
  int nMax;
  /// Maximum allowed value for the quantum number l
  int lMax;
  /// Maximum allowed value for the quantum number m
  int mMax;
  // Operations
public:
  void evalWaveFunctions1D (arma::vec r);
  void evalWaveFunctions3D (arma::mat r);
};

