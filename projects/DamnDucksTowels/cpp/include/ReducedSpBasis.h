#ifndef REDUCEDSPBASIS_H
#define REDUCEDSPBASIS_H

#include "SpBasis.h"

/// class ReducedSpBasis - Implements a single-particle harmonic oscillator basis only considering S-wave
class ReducedSpBasis : public SpBasis
{
  // Operations
public:
  /// Constructor of the ReducedSpBasis class
  ReducedSpBasis (double _omega, int _nMax) : SpBasis(_omega, _nMax, 0)
  {
    type = "ReducedSpBasis";
  };
  /// Destructor
  ~ReducedSpBasis ()
  {
  };
};

#endif
