#include "HartreeFock.h"
#include "quadrature.h"
#include "SpBasis.h"
#include "global.h"

HartreeFock::HartreeFock(System & _system) : Solver(_system, (unsigned int)1), D(_system.particleNumbers.n_rows), occ(_system.particleNumbers.n_rows)
{
  int basisSize = _system.basis->size;
  Basis &basis = *system->basis;
  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  {
    D(pType).eye(basisSize, basisSize);
    occ(pType) = arma::zeros<arma::vec>(basisSize);
    for(unsigned int state = 0; state < _system.particleNumbers(pType); state++)
    {
      occ(pType)(state) = 1;
    }
  }
  for (int bra = 0; bra < basis.size; bra++)
  for (int ket = 0; ket < basis.size; ket++)
  for (int bra2 = 0; bra2 < basis.size; bra2++)
  for (int ket2 = 0; ket2 < basis.size; ket2++)
  {
//    V(bra, ket)(bra2, ket2) = system->inter->get()
  }
}

HartreeFock::~HartreeFock()
{
}

void HartreeFock::run()
{
  int nb_state = system->basis->size;
  arma::field<arma::mat> &H = system->H;
  arma::field<arma::mat> &R = system->R;
  for (unsigned int pType = 0; pType < system->particleNumbers.n_rows; pType++)
  {
    // Temporary vectors and matrices to store eigenvecs and energies
    arma::vec old_indivE = indivEnergies.row(pType).t();
    arma::vec new_indivE;

    // Hamiltonian diagonalization to extract D and e
    arma::eig_sym(new_indivE, D(pType), H(0, pType));

    // Extraction of eigenenergies' sequence
    arma::uvec sorted_ind = arma::sort_index(new_indivE);
    sorted_ind = sorted_ind.head(system->particleNumbers(pType));
    occ(pType).zeros();
    occ(pType)(sorted_ind) = arma::ones<arma::vec>(sorted_ind.n_rows);

    // New derivation of rho
    R(0,pType) = D(pType).t() * arma::diagmat(occ(pType)) * D(pType);

    indivEnergies.row(pType) = new_indivE.t();

    // Convergence check
    cvg = arma::accu(arma::abs(new_indivE - old_indivE)) / nb_state;
  }
}

void HartreeFock::calcH()
{
  Basis &basis = *system->basis;
  arma::field<arma::mat> &H = system->H;
  Interaction &inter = *system->inter;
  int pNum = system->particleNumbers.n_elem;
  if(basis.type == "SpBasis" || basis.type == "ReducedSpBasis")
  {
    SpBasis &spBasis = static_cast<SpBasis&>(basis);
    
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
          potential = 0.5 * NUCLEON_MASS * spBasis.omega * spBasis.omega * arma::accu(wi_p % fun);
          
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
}
