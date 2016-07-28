#include "Solver.h"
#include <sstream>
#include <string>


Solver::Solver(System & _system) : system(&_system),
                indivEnergies(_system.particleNumbers.n_rows, _system.basis->size, arma::fill::zeros), cvg(1000000)
{
}

Solver::~Solver(){}

void Solver::initH (int type)
{
  (void) type;
  int pNum = system->particleNumbers.n_elem;

  for(int pType = 0; pType < pNum; pType++)
  {
    system->Kinetic(pType).eye();
  }
}



std::string Solver::info()
{
  std::stringstream info;
  arma::field<arma::mat> &R = system->R;
  info << "[" << system->name << "] ";
  for (unsigned int pType = 0; pType < R.n_cols; pType++)
  {
    info << "tr(" << system->particleNames[pType] << "): ( ";
    for (unsigned int dType = 0; dType < R.n_cols; dType++)
    {
      info << arma::trace(R(dType, pType)) << " ";
    }
    info << " ), ";
  }
  info << "cvg: " << cvg;
  return info.str();
}

std::string Solver::toString()
{
  return std::string();
}
