#include "MinnesotaS0.h"
#include "quadrature.h"

#include <iostream>

MinnesotaS0::MinnesotaS0(ReducedSpBasis &_basis) : Interaction(_basis, 1), TBME(_basis.size,_basis.size), basis(_basis)
{ 
  V0R = +200.0;
  V0t = -178.0;
  V0s = -91.85;
  
  kR = 1.487;
  kt = 0.639;
  ks = 0.465;
  
  basis = _basis;
  for (int i = 0; i < _basis.size; i++)
    for (int j = 0; j < _basis.size; j++)
      TBME(i,j) = arma::zeros(_basis.size,_basis.size);
    
  calculateTBME(150);
  printTBME();
}

MinnesotaS0::~MinnesotaS0()
{
}

double MinnesotaS0::get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId)
{
  return 0.0;
}

void MinnesotaS0::calculateTBME(int quadratureOrder)
{
  arma::vec wi_p;
  arma::vec pi_p;
  arma::vec fun;
  GET_LAG_ROOTS(quadratureOrder,pi_p,wi_p);
  
  arma::mat wfMatrix;
  basis.evalRadialWaveFunction(wfMatrix, pi_p);
  
  for (int idx1 = 0; idx1 < basis.size; idx1++)
    for (int idx2 = 0; idx2 < basis.size; idx2++)
      for (int idx3 = 0; idx3 < basis.size; idx3++)
	for (int idx4 = 0; idx4 < basis.size; idx4++)
	{
	  double integral = 0.0;
	  int deltas = basis.deltaSpin(idx1,idx3)*basis.deltaSpin(idx2,idx4)-basis.deltaSpin(idx1,idx4)*basis.deltaSpin(idx2,idx3);
	  if (deltas != 0)
	    for (unsigned int r1 = 0; r1 < pi_p.n_elem; r1++)
	    {
	      fun = 0.25 * V0R / kR * pi_p(r1) * (wfMatrix.col(idx1))(r1) * (wfMatrix.col(idx3))(r1) * exp(pi_p(r1)) * wi_p % pi_p  % wfMatrix.col(idx2) % wfMatrix.col(idx4) % ( arma::exp(-kR*arma::pow(pi_p - pi_p(r1),2)) - arma::exp(-kR*arma::pow(pi_p + pi_p(r1),2)) ) % arma::exp(pi_p);
	      integral += 0.5 * deltas * wi_p(r1) * arma::accu(fun);
	      
	      fun = 0.25 * V0s / ks * pi_p(r1) * (wfMatrix.col(idx1))(r1) * (wfMatrix.col(idx3))(r1) * exp(pi_p(r1)) * wi_p % pi_p  % wfMatrix.col(idx2) % wfMatrix.col(idx4) % ( arma::exp(-ks*arma::pow(pi_p - pi_p(r1),2)) - arma::exp(-ks*arma::pow(pi_p + pi_p(r1),2)) ) % arma::exp(pi_p);
	      integral += 0.5 * deltas * wi_p(r1) * arma::accu(fun);
	      
	      fun = 0.25 * V0R / kR * pi_p(r1) * (wfMatrix.col(idx1))(r1) * (wfMatrix.col(idx4))(r1) * exp(pi_p(r1)) * wi_p % pi_p  % wfMatrix.col(idx2) % wfMatrix.col(idx3) % ( arma::exp(-kR*arma::pow(pi_p - pi_p(r1),2)) - arma::exp(-kR*arma::pow(pi_p + pi_p(r1),2)) ) % arma::exp(pi_p);
	      integral += 0.5 * deltas * wi_p(r1) * arma::accu(fun);
	      
	      fun = 0.25 * V0s / ks * pi_p(r1) * (wfMatrix.col(idx1))(r1) * (wfMatrix.col(idx4))(r1) * exp(pi_p(r1)) * wi_p % pi_p  % wfMatrix.col(idx2) % wfMatrix.col(idx3) % ( arma::exp(-ks*arma::pow(pi_p - pi_p(r1),2)) - arma::exp(-ks*arma::pow(pi_p + pi_p(r1),2)) ) % arma::exp(pi_p);
	      integral += 0.5 * deltas * wi_p(r1) * arma::accu(fun);
	    }
	  TBME(idx1,idx2)(idx3,idx4) = integral;
	}
}

void MinnesotaS0::printTBME()
{
  //Printing the names of quantum numbers
  for (int j = 0; j < basis.qNumSize; j++)
      std::cout << basis.qNames[j] << "\t";
  std::cout << "TBME\n";
  for (int idx1 = 0; idx1 < basis.size; idx1++)
    for (int idx2 = 0; idx2 < basis.size; idx2++)
      for (int idx3 = 0; idx3 < basis.size; idx3++)
	for (int idx4 = 0; idx4 < basis.size; idx4++)
	{
	  //Printing quantum number for each basis state
	  std::cout << basis.qNumbers(idx1,0) << "\t";
	  std::cout << basis.qNumbers(idx2,0) << "\t";
	  std::cout << basis.qNumbers(idx3,0) << "\t";
	  std::cout << basis.qNumbers(idx4,0) << "\t";
	  std::cout << TBME(idx1,idx2)(idx3,idx4) << "\n";
	}
}
