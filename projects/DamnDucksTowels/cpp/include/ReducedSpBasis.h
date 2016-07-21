#ifndef REDUCEDSPBASIS_H
#define REDUCEDSPBASIS_H

#include "SpBasis.h"

/// class ReducedSpBasis - 
class ReducedSpBasis : public SpBasis {
  // Attributes
public:
  double omega;
  int lMax;
  int mMax;
  // Operations
public:
  ~ReducedSpBasis ();
};

#endif
