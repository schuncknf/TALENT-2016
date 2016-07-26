#include "HartreeFockBogo.h"
#include <armadillo>

HartreeFockBogo::HartreeFockBogo(System & _system) : Solver(_system, (unsigned int)2), D(_system.particleNumbers.n_rows) 
{
}

HartreeFockBogo::~HartreeFockBogo()
{
}

void HartreeFockBogo::run()
{
}

void HartreeFockBogo::calcH()
{

}
