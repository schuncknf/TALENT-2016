#include "Solver.h"

Solver::Solver(int matSize) : indivEnergies(matSize, arma::fill::zeros), cvg(0)
{
}

