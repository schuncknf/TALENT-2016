#include "RawInteraction.h"

RawInteraction::RawInteraction(Basis & _basis, int _nParticleTypes) :
  Interaction( _basis, _nParticleTypes ),
  potential(_nParticleTypes, _basis.size),
  nParticleTypes(_nParticleTypes)
{
  for (int pType0 = 0; pType0 < nParticleTypes; pType0++)
  for (int state0 = 0; state0 < _basis.size;    state0++)
  {
    potential(pType0, state0) = arma::field<arma::field<arma::mat> >(_nParticleTypes, _basis.size);
    for (int pType1 = 0; pType1 < nParticleTypes; pType1++)
    for (int state1 = 0; state1 < _basis.size;    state1++)
    {
      potential(pType0, state0)(pType1, state1) = arma::field<arma::mat>(_nParticleTypes, _basis.size);
    }
  }
}

RawInteraction::~RawInteraction()
{
}

double & RawInteraction::set(int ti0, int i0, int ti1, int i1, int tj0, int j0, int tj1, int j1)
{
  if(potential(ti0, i0)(ti1, i1)(tj0, j0).empty())
  {
    potential(ti0, i0)(ti1, i1)(tj0, j0) = arma::zeros(nParticleTypes, particleBases(0)->size);
  }
  return potential(ti0, i0)(ti1, i1)(tj0, j0)(tj1, j1);
}

double RawInteraction::get(arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId)
{
  (void) R;
  if(potential(bType(0), bId(0))(bType(1), bId(1))(kType(0), kId(0)).empty())
  {
    return 0.;
  }
  return potential(bType(0), bId(0))(bType(1), bId(1))(kType(0), kId(0))(kType(1), kId(1));
}
