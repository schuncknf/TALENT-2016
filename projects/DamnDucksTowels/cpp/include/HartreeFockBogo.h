#ifndef HARTREEFOCKBOGO_h
#define HARTREEFOCKBOGO_h

#include "Solver.h"
#include "System.h"
#include <armadillo>

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
  void calc (arma::field<arma::mat> & H);
  ~HartreeFockBogo ();
};

#endif
