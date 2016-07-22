#include "HartreeFock.h"

HartreeFock::HartreeFock(System & system) : Solver(system.basis->qNumbers.n_rows), D(system.particleNumbers(0))
	{
	D(0).eye(system.basis->qNumbers.n_rows, system.basis->qNumbers.n_rows);
	RG(1, system.particleNumbers.n_elem);
	RG(0,0).zeros(system.basis->qNumbers.n_rows, system.basis->qNumbers.n_rows);
	for(int i=0; i< system.particleNumbers(0); i++){
		RG(0,0)(i,i,1);
	}
}

void HartreeFock::calc(arma::field<arma::mat> H) {

	int nb_state = system->basis->qNumbers.n_rows;
	// Hamiltonian diagonalization to extract D and e
	arma::eig_sym(indivEnergies,D(0),H(0));

	// Extraction of eigenenergies' sequence
	arma::uvec sorted_ind = sort_index(indivEnergies);

	// Temporary vectors and matrices to store eigenvecs and energies
    arma::vec old_indivE(nb_state);
    old_indivE = indivEnergies ;
	arma::vec new_indivE(nb_state);
	arma::field<arma::mat> new_D(system->particleNumbers(0));
	new_indivE.zeros();
	new_D(0).zeros(nb_state);

	// Eigenstates reorganized with growing energy
	for (int i=0; i < nb_state; i++) {
		int ind = sorted_ind(i);
		new_indivE(i) = indivEnergies(ind);
		new_D(0).col(i) = D(0).col(ind);
	}

	// Eigenvecs and eigenenergies stored in class' attributes
	indivEnergies = new_indivE;
	D(0) = new_D(0);

	// New derivation of rho
	RG(0,0) = D(0)*D(0);

    // Convergence check
    cvg = 0;
    for(int i=0; i < nb_state; i++) {
        cvg += abs(indivEnergies(i) - old_indivE(i));
    }
    cvg = cvg / nb_state;
}

