#include "HartreeFockBogo.h"
#include <armadillo>

HartreeFockBogo::HartreeFockBogo(System &_system) : Solver(_system), D(_system.particleNumbers.n_rows)
{
  int basisSize = _system.basis->size;

  for (unsigned int pType = 0; pType < _system.particleNumbers.n_rows; pType++)
  {
    D(pType).eye(2 * basisSize, 2 * basisSize);
    occ(pType) = arma::zeros<arma::vec>(2 * basisSize);

    for (unsigned int state = 0; state < _system.particleNumbers(pType); state++)
    {
      occ(pType)(state) = 1;
    }
  }
}

HartreeFockBogo::~HartreeFockBogo()
{
}

void HartreeFockBogo::run()
{
  int nb_state = system->basis->size;
  arma::field<arma::mat> &R = system->R;

  for (unsigned int pType = 0; pType < system->particleNumbers.n_rows; pType++)
  {
    arma::mat H = system->Kinetic(pType) + system->Gamma(pType);
    // Temporary vectors and matrices to store eigenvecs and energies
    arma::vec old_indivE = indivEnergies.row(pType).t();
    arma::vec new_indivE;
    // "Generalized" hamiltonian matrix
    arma::mat bigH( 2 * nb_state, 2 * nb_state );
    bigH.submat( 0, 0, nb_state - 1, nb_state - 1 ) = H ; //bigH(0,0)
    bigH.submat( 0, nb_state, nb_state - 1, 2 * nb_state - 1 ).zeros() ; //bigH(0,1)
    bigH.submat( nb_state, 0, 2 * nb_state - 1, nb_state - 1 ).zeros() ; //bigH(1,0)
    bigH.submat( nb_state, nb_state, 2 * nb_state - 1, 2 * nb_state - 1 ) = - H; //BigH(1,1)
    // Hamiltonian diagonalization to extract D and e
    arma::eig_sym(new_indivE, D(pType), bigH);
    // Extraction of eigenenergies' sequence
    arma::uvec sorted_ind = arma::sort_index(new_indivE);
    sorted_ind = sorted_ind.head(system->particleNumbers(pType));
    occ(pType).zeros();
    occ(pType)(sorted_ind) = arma::ones<arma::vec>(sorted_ind.n_rows);
    // New derivation of rho and kappa
    arma::mat U = D(pType).submat(0, 0, nb_state - 1, nb_state - 1).t();
    arma::mat V = D(pType).submat(0, nb_state, nb_state - 1, 2 * nb_state - 1).t();
    R(0, pType) = V * V.t();
    R(1, pType) = V * U.t();
    indivEnergies.row(pType) = new_indivE.t();
    // Convergence check
    cvg = arma::accu(arma::abs(new_indivE - old_indivE)) / nb_state;
  }
}

void HartreeFockBogo::calcH()
{
}

