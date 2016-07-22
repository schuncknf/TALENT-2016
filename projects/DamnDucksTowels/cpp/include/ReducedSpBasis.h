#ifndef REDUCEDSPBASIS_H
#define REDUCEDSPBASIS_H

#include "SpBasis.h"

/// class ReducedSpBasis 
class ReducedSpBasis : public SpBasis {
public:
  ReducedSpBasis(double _omega, int _nMax) : SpBasis(_omega, _nMax, 0)
  {
    type = "ReducedSpBasis";
  };
  ~ReducedSpBasis ()
  {
  };
};

#endif
