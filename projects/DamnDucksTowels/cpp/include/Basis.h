#ifndef BASIS_H
#define BASIS_H

#include <armadillo>
#include <vector>
#include <string>

/// class Basis - 
class Basis {
  // Attributes
public:
  /// Type of basis (S-waves only or not, etc)
  std::string type;
  /// Number of quantum numbers for each states
  int qNumSize;
  /// The k-th element of this vector is the name of the k-th quantum number (e.g. "n
  std::vector<std::string> qNames;
  /// The k-th column of the i-th row of this matrix contain the k-th quantum number 
  arma::imat qNumbers;
  /// Size of the basis
  int basisSize;
  // Operations
  /// Printing basis states with quantum numbers
  void printBasis();
public:
  Basis (std::string type, std::vector<std::string> qNames);
  void calcOverlaps (arma::mat & overlaps);
  virtual ~Basis () = 0;
};

#endif
