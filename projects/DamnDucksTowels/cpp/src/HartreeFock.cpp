#include "HartreeFock.h"

HartreeFock::HartreeFock(System & system) 
	: D(system.particleNumbers(0),nb_state,arma::fill::eye), 
	rho(nb_state,nb_state,arma::fill::zeros),
	e(nb_state,arma::fill::zeros),
	
	{
	for(int i=0; i< system.particleNumbers(0); i++){
		rho(i,i,1);
	}
	RG(0,0) = rho;
	}

void HartreeFock::iter(arma::field<arma::mat> H) {
	// Hamiltonian diagonalization to extract D and e
	arma::eig_gen(e,D,H(0));

	// Extraction of eigenenergies' sequence
	arma::uvec sorted_ind = sort_index(e);

	// Temporary vectors and matrices to store eigenvecs and energies
	arma::vec new_e(nb_state,arma::fill::zeros);
	arma::mat new_D(nb_state,nb_state,arma::fill::zeros);

	// Eigenstates reorganized with growing energy
	for (int i=0; i < nb_state; i++) {
		int ind = sorted_ind(i);
		new_e(i) = e(ind);
		new_D.col(i) = D.col(ind);
	}

	// Eigenvecs and eigenenergies stored in class' attributes
	e = new_e;
	D = new_D;

	// New derivation of rho
	rho = D*D;
}

