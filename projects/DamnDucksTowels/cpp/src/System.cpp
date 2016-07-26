#include "System.h"

System::System (std::string _name, Basis & _basis, arma::ivec _particleNumbers, std::vector<std::string> _particleNames, Interaction & _inter) :
    name(_name), basis(&_basis), particleNumbers(_particleNumbers), particleNames(_particleNames), inter(&_inter) {}

System::~System(){}

arma::mat & System::getH(int dType, int pType)
{
  return H(dType, pType);
}

arma::mat & System::getR(int dType, int pType)
{
  return R(dType, pType);
}

std::string System::info()
{
  return std::string();
}

std::string System::toString()
{
  return std::string();
}

void System::addParticleType(std::string name, int number)
{
  particleNumbers.resize(particleNumbers.n_rows + 1);
  particleNumbers(particleNumbers.n_rows - 1) = number;
  particleNames.push_back(name);
}
