#include "Nucleus.h"

Nucleus::Nucleus (Basis & basis, Interaction & _inter, int nbNeut, int nbProt, int _nPoints) 
    : System("Nucleus", basis, arma::ivec {nbNeut,nbProt}, std::vector<std::string> {"Neutrons","Protons"}, _inter, _nPoints)
{
}

Nucleus::~Nucleus()
{
}

void Neutron::calcH()
{

}
