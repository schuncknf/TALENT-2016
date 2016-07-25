#include "HartreeFock.h"

HartreeFock::HartreeFock(System & _system) : Solver(_system, (unsigned int)1), D(_system.particleNumbers.n_rows), occ(_system.particleNumbers.n_rows)
{
  int basisSize = _system.basis->size;
  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  {
    D(pType).eye(basisSize, basisSize);
    occ(pType) = arma::zeros<arma::vec>(basisSize);
    for(unsigned int state = 0; state < _system.particleNumbers(pType); state++)
    {
      occ(pType)(state) = 1;
    }
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
  
}
