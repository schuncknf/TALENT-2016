#include <armadillo>

#include "Interaction.h"

Interaction::Interaction(arma::field<Basis *> _particleBases) : particleBases(_particleBases)
{
}

Interaction::Interaction(Basis & _basis, int _nParticleTypes)
{
  particleBases = arma::field<Basis *>(_nParticleTypes);
  for (int pType = 0; pType < _nParticleTypes; pType++)
  {
    particleBases(pType) = &_basis;
  }
}

Interaction::~Interaction()
{
}

double Interaction::get(arma::field<arma::mat> & R, int ti, int i, int tj, int j)
{
  arma::ivec bType ( { ti } );
  arma::ivec bId   ( { i  } );
  arma::ivec kType ( { tj } );
  arma::ivec kId   ( { j  } );
  return get(R, bType, bId, kType, kId);
}

double Interaction::get(arma::field<arma::mat> & R, int ti0, int i0, int ti1, int i1, int tj0, int j0, int tj1, int j1)
{
  arma::ivec bType ( { ti0, ti1 } );
  arma::ivec bId   ( { i0 , i1  } );
  arma::ivec kType ( { tj0, tj1 } );
  arma::ivec kId   ( { j0 , j1  } );
  return get(R, bType, bId, kType, kId);
}
