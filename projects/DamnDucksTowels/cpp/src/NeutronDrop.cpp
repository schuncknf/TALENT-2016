#include <armadillo>
#include <string>
#include <vector>

#include "quadrature.h"
#include "global.h"
#include "NeutronDrop.h"
#include "SpBasis.h"
#include "FullSpBasis.h"


NeutronDrop::NeutronDrop(Basis &_basis, Interaction &_inter, int _nbNeut, double _omega, int _nPoints) :
  System(std::string("NeutronDrop"), _basis, arma::ivec(
{
  _nbNeut
}), std::vector<std::string>({"neutron"}), _inter, _nPoints), omega(_omega),TBME(_basis.size, _basis.size)

{
  if(basis->type == "SpBasis" || basis->type == "ReducedSpBasis")
  {
     for (int i = 0; i < basis->size; i++)
      for (int j = 0; j < basis->size; j++)
      {
        TBME(i, j) = arma::zeros(_basis.size, basis->size);
        for (int k = 0; k < basis->size; k++)
          for (int l = 0; l < basis->size; l++)
          {
            arma::field<arma::mat> dummyR;
            TBME(i, j)(k, l) = inter->get(dummyR, 0, i, 0, l, 0, j, 0, k);
          }
      }
  }
  else if (basis->type == "FullSpBasis")
  {
    for (int i = 0; i < basis->size; i++)
      for (int j = 0; j < basis->size; j++)
      {
        if (_basis.qNumbers(i, 1) != _basis.qNumbers(j, 1)) continue;
        if (_basis.qNumbers(i, 2) != _basis.qNumbers(j, 2)) continue;
        TBME(i, j) = arma::zeros(_basis.size, basis->size);
        
        for (int k = 0; k < basis->size; k++)
          for (int l = 0; l < basis->size; l++)
          {
            if (_basis.qNumbers(k, 1) != _basis.qNumbers(l, 1)) continue;
            if (_basis.qNumbers(k, 2) != _basis.qNumbers(l, 2)) continue;
            arma::field<arma::mat> dummyR;
            TBME(i, j)(k, l) = inter->get(dummyR, 0, i, 0, l, 0, j, 0, k);
          }
      }
  }
}

NeutronDrop::~NeutronDrop()
{
}

void NeutronDrop::calcH()
{
  if (basis->type == "SpBasis" || basis->type == "ReducedSpBasis")
  {
    SpBasis &spBasis = static_cast<SpBasis &>(*basis);
    //Initialization of quadrature
    arma::vec wi_p;
    arma::vec pi_p;
    arma::vec fun;
    GET_LAG_ROOTS(nPoints, pi_p, wi_p);
    arma::mat wfMatrix;
    arma::mat wfDerMatrix;
    spBasis.evalRadialWaveFunction(wfMatrix, pi_p);
    spBasis.evalDerivativeRadialWaveFunction(wfDerMatrix, pi_p);

    for (int bra = 0; bra < spBasis.size; bra++)
      for (int ket = 0; ket < spBasis.size; ket++)
      {
        int lBra = spBasis.qNumbers(bra, 1);
        int lKet = spBasis.qNumbers(ket, 1);
        int mBra = spBasis.qNumbers(bra, 2);
        int mKet = spBasis.qNumbers(ket, 2);
        int sBra = spBasis.qNumbers(bra, 3);
        int sKet = spBasis.qNumbers(ket, 3);

        if ((lBra == lKet) && (mBra == mKet) && (sBra == sKet))
        {
          // Calculation of the kinetic part
          double kinetic, potential;
          fun = ( arma::pow(pi_p, 2) % wfDerMatrix.col(bra) % wfDerMatrix.col(ket) + lBra * (lBra + 1) * wfMatrix.col(bra) % wfMatrix.col(ket)) % arma::exp(pi_p);
          kinetic = HBAR2_2M * arma::accu(wi_p % fun);
          //kinetic = HBAR*HBAR/2.0/NUCLEON_MASS * arma::accu(wi_p % fun);
          // Calculation of the harmonic part
          fun = arma::pow(pi_p, 4) % wfMatrix.col(bra) % wfMatrix.col(ket) % arma::exp(pi_p);
          potential = 0.5 * NUCLEON_MASS * pow(spBasis.omega, 2) * arma::accu(wi_p % fun);

          if (abs(omega - spBasis.omega) > 1.0e-7)
          {
            Kinetic(0)(bra, ket) = kinetic + potential;
          }
          else
          {
            Kinetic(0)(bra, ket) = 0;

            if (bra == ket)
            {
              Kinetic(0)(bra, ket) = 10 * (2 * basis->qNumbers(bra, 0) + 3. / 2.);
            }
          }
        }
        else
        {
          Kinetic(0)(bra, ket) = 0.0;
        }
      }
  }
  else if (basis->type == "FullSpBasis")
  {
    FullSpBasis &spBasis = static_cast<FullSpBasis &>(*basis);
    //Initialization of quadrature
    arma::vec wi_p;
    arma::vec pi_p;
    arma::vec fun;
    GET_LAG_ROOTS(nPoints, pi_p, wi_p);
    arma::mat wfMatrix;
    arma::mat wfDerMatrix;
    spBasis.evalRadialWaveFunction(wfMatrix, pi_p);
    spBasis.evalDerivativeRadialWaveFunction(wfDerMatrix, pi_p);

    for (int bra = 0; bra < spBasis.size; bra++)
      for (int ket = 0; ket < spBasis.size; ket++)
      {
        int lBra = spBasis.qNumbers(bra, 1);
        int lKet = spBasis.qNumbers(ket, 1);
        int jBra = spBasis.qNumbers(bra, 2);
        int jKet = spBasis.qNumbers(ket, 2);

        if ((lBra == lKet) && (jBra == jKet))
        {
          if (abs(omega - spBasis.omega) > 1.0e-7)
          {
            // Calculation of the kinetic part
            double kinetic, potential;
            fun = ( arma::pow(pi_p, 2) % wfDerMatrix.col(bra) % wfDerMatrix.col(ket) + lBra * (lBra + 1) * wfMatrix.col(bra) % wfMatrix.col(ket)) % arma::exp(pi_p);
            kinetic = HBAR2_2M * arma::accu(wi_p % fun);
            //kinetic = HBAR*HBAR/2.0/NUCLEON_MASS * arma::accu(wi_p % fun);
            // Calculation of the harmonic part
            fun = arma::pow(pi_p, 4) % wfMatrix.col(bra) % wfMatrix.col(ket) % arma::exp(pi_p);
            potential = 0.5 * NUCLEON_MASS * pow(spBasis.omega, 2) * arma::accu(wi_p % fun);

            Kinetic(0)(bra, ket) = kinetic + potential;
          }
          else
          {
            if (bra == ket)
            {
              Kinetic(0)(bra, ket) = HBAR * omega * (2 * basis->qNumbers(bra, 0) + basis->qNumbers(bra, 1) + 3. / 2.);
            }
          }
        }
        else
        {
          Kinetic(0)(bra, ket) = 0.0;
        }
      }
   
    for (int i = 0; i < basis->size; i++)
      for (int j = 0; j < basis->size; j++)
      {
        if (spBasis.qNumbers(i, 1) != spBasis.qNumbers(j, 1)) continue;
        if (spBasis.qNumbers(i, 2) != spBasis.qNumbers(j, 2)) continue;
        Gamma(0)(i, j) = arma::accu(TBME(i, j) % R(0, 0));
      }
  }
  else
  {
    throw std::runtime_error("Unknown basis.");
  }

}
