#include "HartreeFockBogo.h"

HartreeFockBogo::HartreeFockBogo(System & system) : Solver(system.basis->qNumbers.n_rows), D(system.particleNumbers(0)) 
{
	D(0).eye(system.basis->qNumbers.n_rows, system.basis->qNumbers.n_rows);
	RG(2, system.particleNumbers.n_elem);
	RG(0,0).zeros(system.basis->qNumbers.n_rows, system.basis->qNumbers.n_rows);
	RG(1,0).zeros(system.basis->qNumbers.n_rows, system.basis->qNumbers.n_rows);
	for(int i=0; i < system.particleNumbers(0); i++) {
		RG(0, 0)(i, i) = 1;
	}
}

void HartreeFockBogo::calc(arma::field<arma::mat> & H) {
  (void) H;
}
