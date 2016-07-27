#include <armadillo>
#include <string>
#include <vector>

#include "quadrature.h"
#include "global.h"
#include "NeutronDrop.h"
#include "SpBasis.h"


NeutronDrop::NeutronDrop(Basis & _basis, Interaction & _inter, int _nbNeut, int _nPoints) :
  System(std::string("NeutronDrop"), _basis, arma::ivec({_nbNeut}), std::vector<std::string>({"neutron"}), _inter, _nPoints)
{
}

NeutronDrop::~NeutronDrop()
{
}

void NeutronDrop::calcH()
{
  int pNum = particleNumbers.n_elem;
  if(basis->type == "SpBasis" || basis->type == "ReducedSpBasis")
  {
    SpBasis &spBasis = static_cast<SpBasis&>(*basis);
    
    //Initialization of quadrature
    arma::vec wi_p;
    arma::vec pi_p;
    arma::vec fun;
    int nodesNum = 90;
    GET_LAG_ROOTS(nodesNum,pi_p,wi_p);
    arma::mat wfMatrix;
    arma::mat wfDerMatrix;
    spBasis.evalRadialWaveFunction(wfMatrix, pi_p);
    spBasis.evalDerivativeRadialWaveFunction(wfDerMatrix, pi_p);
    for (int pTyp = 0; pTyp < pNum; pTyp++)
    {
      arma::mat gamma(spBasis.size, spBasis.size, arma::fill::zeros);
      for (int bra = 0; bra < spBasis.size; bra++)
      for (int ket = 0; ket < spBasis.size; ket++)
      {
        int lBra = spBasis.qNumbers(bra,1);
        int lKet = spBasis.qNumbers(ket,1);
        int mBra = spBasis.qNumbers(bra,2);
        int mKet = spBasis.qNumbers(ket,2);
        int sBra = spBasis.qNumbers(bra,3);
        int sKet = spBasis.qNumbers(ket,3);
        if ((lBra == lKet) && (mBra == mKet) && (sBra == sKet))
        {
          // Calculation of the kinetic part
          double kinetic, potential;
          fun = ( arma::pow(pi_p,2) % wfDerMatrix.col(bra) % wfDerMatrix.col(ket) + lBra*(lBra+1) * wfMatrix.col(bra) % wfMatrix.col(ket)) % arma::exp(pi_p);
          kinetic = HBAR*HBAR/2.0/NUCLEON_MASS * arma::accu(wi_p % fun);
        
            // Calculation of the harmonic part
            fun = arma::pow(pi_p,4) % wfMatrix.col(bra) % wfMatrix.col(ket) % arma::exp(pi_p);
            potential = 0.5 * NUCLEON_MASS * pow(spBasis.omega,2) * arma::accu(wi_p % fun);
            
            H(0,pTyp)(bra,ket) = kinetic + potential;
        }
        else
        {
          H(0,pTyp)(bra,ket) = 0.0;
        }
      }
    }
  }
  else
  {
    throw std::runtime_error("Unknown basis.");
  }
  for (int pTyp = 0; pTyp < pNum; pTyp++)
  {
    arma::field<arma::mat> TBME = arma::field<arma::mat>(basis->size, basis->size);
    for (int i = 0; i < basis->size; i++)
    for (int j = 0; j < basis->size; j++)
    {
      TBME(i,j) = arma::zeros(basis->size,basis->size);
      for (int k = 0; k < basis->size; k++)
      for (int l = 0; l < basis->size; l++)
      {
        arma::field<arma::mat> dummyR;
        TBME(i,j)(k,l) = inter->get(dummyR, pTyp, i, pTyp, l, pTyp, j, pTyp, k);
      }
    }
    TBME(0,0).print();
    for (int bra = 0; bra < basis->size; bra++)
    for (int ket = 0; ket < basis->size; ket++)
    {
      H(0,pTyp)(bra,ket) += arma::accu(TBME(bra,ket) % R(0,pTyp));
    }
  }


}
