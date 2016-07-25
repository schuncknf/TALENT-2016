#ifndef BASIS_H
#define BASIS_H

#include <armadillo>
#include <vector>
#include <string>

/// class Basis - Abstract class representing a basis
class Basis {
  // Attributes
public:
  /// Type of basis (S-waves only or not, etc)
  std::string type;
  /// The k-th element of this vector is the name of the k-th quantum number (e.g. "n
  std::vector<std::string> qNames;
  /// Number of quantum numbers for each states
  int qNumSize;
  /// The k-th column of the i-th row of this matrix contain the k-th quantum number 
  arma::imat qNumbers;
  /// Size of the basis
  int size;
  // Operations
public:
  /// Constructor
  Basis (std::string _type, std::vector<std::string> _qNames);
  /// Pure virtual destructor
  virtual ~Basis () = 0;
  /// Return few information about the basis
  std::string info ();
  /// Return a lot of information about the basis
  std::string toString ();
};

#endif
