/// class Basis - 
class Basis {
  // Attributes
public:
  /// Type of basis (S-waves only or not)
  int type;
  /// Number of quantum number for each states
  int qNumSize;
  /// The k-th element of this vector is the name of the k-th quantum number (e.g. "n
  std::vector<std::string> qNames;
  /// The k-th column of the i-th row of this matrix contain the k-th quantum number 
  arma::imat qNumbers;
  // Operations
public:
  virtual void calcOverlaps () = 0;
  virtual void calcPOverlaps () = 0;
};

