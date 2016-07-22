#include "Solver.h"

Solver::Solver(int matSize) : indivEnergies(matSize), cvg{0}
{
	indivEnergies.zeros();
}

