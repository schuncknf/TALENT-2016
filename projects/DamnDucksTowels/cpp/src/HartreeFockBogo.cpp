#include "HartreeFockBogo.h"

HartreeFockBogo::HartreeFockBogo(System & system) {
	D(nb_state,nb_state,fill::eye);
	rho(nb_state,nb_state,fill::zeros);
	kappa(nb_state,nb_state,fill:zeros);
	e(nb_state,fill::zeros);
}

void HartreeFockBogo::iter(arma::field<arma::mat> H) {

}
