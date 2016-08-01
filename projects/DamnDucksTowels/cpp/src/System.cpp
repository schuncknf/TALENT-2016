#include "System.h"

System::System (std::string _name, Basis &_basis, arma::ivec _particleNumbers, std::vector<std::string> _particleNames, Interaction &_inter, int _nPoints) :
  name(_name), basis(&_basis), particleNumbers(_particleNumbers), particleNames(_particleNames), inter(&_inter), nPoints(_nPoints)
{
  R = arma::field<arma::mat>(2, particleNumbers.n_rows);
  Kinetic = arma::field<arma::mat>(particleNumbers.n_rows);
  Gamma   = arma::field<arma::mat>(particleNumbers.n_rows);
  Delta   = arma::field<arma::mat>(particleNumbers.n_rows);

  for (unsigned int pType = 0; pType < particleNumbers.n_rows; pType++)
  {
    for (unsigned int dType = 0; dType < 2; dType++)
    {
      R(dType, pType) = arma::zeros(basis->size, basis->size);

      for (unsigned int state = 0; state < particleNumbers(pType); state++)
      {
        R(dType, pType)(state, state) = 1;
      }
    }

    Kinetic(pType) = arma::zeros(basis->size, basis->size);
    Gamma(pType)   = arma::zeros(basis->size, basis->size);
    Delta(pType)   = arma::zeros(basis->size, basis->size);
  }
}

System::~System() {}

arma::mat &System::getKinetic(int pType)
{
  return Kinetic(pType);
}

arma::mat &System::getGamma(int pType)
{
  return Gamma(pType);
}

arma::mat &System::getDelta(int pType)
{
  return Delta(pType);
}

arma::mat &System::getR(int dType, int pType)
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
