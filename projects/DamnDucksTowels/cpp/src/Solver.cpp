#include "Solver.h"
#include <sstream>
#include <string>


Solver::Solver(System & _system, unsigned int _dNumber) : system(&_system),
                indivEnergies(_system.particleNumbers.n_rows, _system.basis->size, arma::fill::zeros), cvg(0)
{
  arma::field<arma::mat> &R = _system.R;
  arma::field<arma::mat> &H = _system.H;
  R = arma::field<arma::mat>(_system.particleNumbers.n_rows, _dNumber);
  H = arma::field<arma::mat>(_system.particleNumbers.n_rows, _dNumber);
  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  for (unsigned int dType = 0; dType < _dNumber; dType++)
  {
    R(pType, dType) = arma::zeros(_system.basis->size, _system.basis->size);
    H(pType, dType) = arma::zeros(_system.basis->size, _system.basis->size);
    for (unsigned int state = 0; state < _system.particleNumbers(pType); state++)
    {
      R(pType, dType) = 1;
    }
  }
}

Solver::~Solver(){}

std::string Solver::info()
{
  std::stringstream info;
  arma::field<arma::mat> &R = system->R;
  info << "[" << system->name << "] ";
  for (unsigned int pType = 0; pType < R.n_rows; pType++)
  {
    info << "tr(" << system->particleNames[pType] << "): ( ";
    for (unsigned int dType = 0; dType < R.n_cols; dType++)
    {
      info << arma::trace(R(pType, dType)) << " ";
    }
    info << " ), ";
  }
  info << "cvg: " << cvg << std::endl;
  return info.str();
}

std::string Solver::toString()
{
  return std::string();
}
