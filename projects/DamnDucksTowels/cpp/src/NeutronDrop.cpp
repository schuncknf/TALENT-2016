#include "NeutronDrop.h"
#include <armadillo>
#include <string>
#include <vector>

NeutronDrop::NeutronDrop(int _nbNeut, Basis & _basis, arma::field<arma::mat> &_TBME) : System(std::string("NeutronDrop"), _basis, arma::ivec({_nbNeut}), std::vector<std::string>({"neutron"}),  _TBME){}
void NeutronDrop::calcH0 (arma::field<arma::mat> & H0, int type)
{}
void NeutronDrop::calcH (arma::field<arma::mat> & H, arma::field<arma::mat> & RG)
{}
void NeutronDrop::calcKineticField(arma::field<arma::mat> & RG)
{}
NeutronDrop::~NeutronDrop()
{}
