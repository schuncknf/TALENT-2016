#include "NeutronDrop.h"
#include <armadillo>
#include <string>
#include <vector>

NeutronDrop::NeutronDrop(Basis & _basis, Interaction & _inter, int _nbNeut, int _nPoints) :
  System(std::string("NeutronDrop"), _basis, arma::ivec({_nbNeut}), std::vector<std::string>({"neutron"}), _inter, _nPoints)
{
}

NeutronDrop::~NeutronDrop()
{
}

void NeutronDrop::calcH()
{
}
