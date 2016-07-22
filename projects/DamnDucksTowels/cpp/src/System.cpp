#include "System.h"

 System::System (std::string _name, Basis & _basis,
                 arma::ivec _particleNumbers, std::vector<std::string> _particleNames,
                 arma::field<arma::mat> & _TBME) : name(_name), basis(&_basis), particleNumbers(_particleNumbers), particleNames(_particleNames),
                 TBME(_TBME) {}
