#include "HartreeFock.h"

HartreeFock::HartreeFock(System & system) 
	{
	int nb_state = system.basis->qNumbers.n_rows;
	D(system.particleNumbers(0), nb_state);
	D.eye();
	RG(2, system.particleNumbers.n_elem);
	indivEnergies(nb_state);
	indivEnergies.zeros();
	for(int i=0; i< system.particleNumbers(0); i++){
		RG(0,0)(i,i,1);
	}
}

void HartreeFock::iter(arma::field<arma::mat> H) {

	int nb_state = system->basis->qNumbers.n_rows;
	// Hamiltonian diagonalization to extract D and e
	arma::eig_sym(indivEnergies,D,H(0));

	// Extraction of eigenenergies' sequence
	arma::uvec sorted_ind = sort_index(indivEnergies);

	// Temporary vectors and matrices to store eigenvecs and energies
	arma::vec new_indivE(nb_state);
	arma::mat new_D(nb_state,nb_state);
	new_indivE.zeros();
	new_D.zeros();

	// Eigenstates reorganized with growing energy
	for (int i=0; i < nb_state; i++) {
		int ind = sorted_ind(i);
		new_indivE(i) = indivEnergies(ind);
		new_D.col(i) = D.col(ind);
	}

	// Eigenvecs and eigenenergies stored in class' attributes
	indivEnergies = new_indivE;
	D = new_D;

	// New derivation of rho
	RG(0,0) = D*D;
}

