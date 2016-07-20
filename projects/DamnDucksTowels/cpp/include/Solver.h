#include "Fields.h"

/// class Solver - 
class Solver {
  // Attributes
public:
  Fields* fields;
  // Operations
public:
  virtual std::string run () = 0;
  virtual std::string info () = 0;
};

