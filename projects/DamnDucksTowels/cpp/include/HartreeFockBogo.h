#ifndef HARTREEFOCKBOGO_h
#define HARTREEFOCKBOGO_h

#include <armadillo>

#include "Solver.h"
#include "System.h"

/// class HartreeFockBogo - 
class HartreeFockBogo : public Solver {
  // Attributes
public:
  /// Transformation matrix from HO to HF
  arma::field<arma::mat> D;
  /// Occupation numbers
  arma::field<arma::vec> occ;
  // Operations
public:
  HartreeFockBogo (System & system);
  ~HartreeFockBogo ();
  void run ();
  void calcH ();
};

#endif
