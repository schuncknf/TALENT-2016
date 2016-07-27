#include "Solver.h"
#include <sstream>
#include <string>


Solver::Solver(System & _system, unsigned int _dNumber) : system(&_system),
                indivEnergies(_system.particleNumbers.n_rows, _system.basis->size, arma::fill::zeros), cvg(1000000)
{
  arma::field<arma::mat> &R = _system.R;
  arma::field<arma::mat> &H = _system.H;
  R = arma::field<arma::mat>(_system.particleNumbers.n_rows, _dNumber);
  H = arma::field<arma::mat>(_system.particleNumbers.n_rows, _dNumber);
  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  for (unsigned int dType = 0; dType < _dNumber; dType++)
  {
    R(dType, pType) = arma::zeros(_system.basis->size, _system.basis->size);
    H(dType, pType) = arma::zeros(_system.basis->size, _system.basis->size);
    for (unsigned int state = 0; state < _system.particleNumbers(pType); state++)
    {
      R(dType, pType) = 1;
    }
  }
}

Solver::~Solver(){}

void Solver::initH (int type)
{
  (void) type;
  int pNum = system->H.n_cols;

  for(int pType = 0; pType < pNum; pType++)
  {
    system->H(0, pType).eye();
  }
}



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
      info << arma::trace(R(dType, pType)) << " ";
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
