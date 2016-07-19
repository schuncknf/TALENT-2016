#include "Fields.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  Fields* fields;
  // Operations
public:
  virtual void run () = 0;
  virtual void info () = 0;
};

