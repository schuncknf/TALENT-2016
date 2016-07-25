#include "NeutronDrop.h"
#include <armadillo>
#include <string>
#include <vector>

NeutronDrop::NeutronDrop(Basis & _basis, Interaction & _inter, int _nbNeut) :
  System(std::string("NeutronDrop"), _basis, arma::ivec({_nbNeut}), std::vector<std::string>({"neutron"}), _inter)
{
}

NeutronDrop::~NeutronDrop()
{
}

void NeutronDrop::calcH0 (int type)
{
  (void) type;
}

void NeutronDrop::calcH ()
{
}
