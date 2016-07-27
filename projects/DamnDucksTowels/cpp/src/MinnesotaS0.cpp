#include "MinnesotaS0.h"
#include "quadrature.h"

#include <sstream>
MinnesotaS0::MinnesotaS0(ReducedSpBasis &_basis, int _nPoints) : Interaction(_basis, 1), TBME(_basis.size,_basis.size), nPoints(_nPoints), basis(_basis)
{   
  basis = _basis;
  for (int i = 0; i < _basis.size; i++)
    for (int j = 0; j < _basis.size; j++)
      TBME(i,j) = arma::zeros(_basis.size,_basis.size);
    
  calc();
  toString();
}

MinnesotaS0::~MinnesotaS0()
{
}

double MinnesotaS0::get (arma::field<arma::mat> & R, arma::ivec & bType, arma::ivec & bId, arma::ivec & kType, arma::ivec & kId)
{
  (void) R;
  if(bType.n_elem != 2 || kType.n_elem != 2 || bId.n_elem != 2 || kId.n_elem != 2)
  {
    return 0.;
  }
  
  return TBME(bId(0), bId(1))(kId(0), kId(1));
}

void MinnesotaS0::calc()
{
  //Parameters of Minnesota potential
  // V0t and kt are not used so they are commented out
  double V0R = +200.0;
  //double V0t = -178.0;
  double V0s = -91.85;
  
  double kR = 1.487;
  //double kt = 0.639;
  double ks = 0.465;
  
  //Gauss-Lagerre quadrature
  arma::vec wi_p;
  arma::vec pi_p;
  GET_LAG_ROOTS(nPoints,pi_p,wi_p);
  
  //Evaluating wave functions
  arma::mat wfMatrix;
  basis.evalRadialWaveFunction(wfMatrix, pi_p);
  
  //Evaluating intermediate results
  arma::field<arma::vec> temp(basis.size,basis.size);
  for (int idx1 = 0; idx1 < basis.size; idx1++)
    for (int idx2 = 0; idx2 < basis.size; idx2++)
      temp(idx1,idx2) = wi_p % pi_p % wfMatrix.col(idx1) % wfMatrix.col(idx2) % exp(pi_p);
  
  arma::mat temp_k = arma::zeros(pi_p.n_elem,pi_p.n_elem);
  for (unsigned int r1 = 0; r1 < pi_p.n_elem; r1++)
    for (unsigned int r2 = 0; r2 < pi_p.n_elem; r2++)
    {
      temp_k(r1,r2) = 0.25 * V0R / kR * (exp(-kR*pow(pi_p(r2) - pi_p(r1),2)) - exp(-kR*pow(pi_p(r2) + pi_p(r1),2))) + 0.25 * V0s / ks * (exp(-ks*pow(pi_p(r2) - pi_p(r1),2)) - exp(-ks*pow(pi_p(r2) + pi_p(r1),2)));
    }
  
  //Calculating anti-symmetrized two-body matrix elements
  arma::mat fun;
  for (int idx1 = 0; idx1 < basis.size; idx1++)
    for (int idx2 = 0; idx2 < idx1; idx2++)
      for (int idx3 = 0; idx3 <= idx1; idx3++)
	for (int idx4 = 0; (idx4 <= idx2) || (idx4 < idx3); idx4++)
	{
	  double integral = 0.0;
	  //Delta of spins
	  int deltas = basis.deltaSpin(idx1,idx3)*basis.deltaSpin(idx2,idx4)-basis.deltaSpin(idx1,idx4)*basis.deltaSpin(idx2,idx3);
	  if (deltas != 0)
	    {
	      fun = (temp(idx1,idx3) * temp(idx2,idx4).t()) % temp_k;
	      integral += 0.5 * deltas * arma::accu(fun);
	      
	      fun = (temp(idx1,idx4) * temp(idx2,idx3).t()) % temp_k;
	      integral += 0.5 * deltas * arma::accu(fun);
	    }
	  //else integral is 0.0 as preset;  

	  TBME(idx1,idx2)(idx3,idx4) = +integral;
	  TBME(idx2,idx1)(idx3,idx4) = -integral;
	  TBME(idx1,idx2)(idx4,idx3) = -integral;
	  TBME(idx2,idx1)(idx4,idx3) = +integral;
	  
	  TBME(idx3,idx4)(idx1,idx2) = +integral;
	  TBME(idx4,idx3)(idx1,idx2) = -integral;
	  TBME(idx3,idx4)(idx2,idx1) = -integral;
	  TBME(idx4,idx3)(idx2,idx1) = +integral;
	}
}

std::string MinnesotaS0::info()
{
  return std::string();
}

std::string MinnesotaS0::toString()
{
  std::stringstream ss;
  //Printing the names of quantum numbers
  for (int j = 0; j < basis.qNumSize; j++)
      ss << basis.qNames[j] << "\t";
  ss << "TBME\n";
  for (int idx1 = 0; idx1 < basis.size; idx1++)
    for (int idx2 = 0; idx2 < basis.size; idx2++)
      for (int idx3 = 0; idx3 < basis.size; idx3++)
	for (int idx4 = 0; idx4 < basis.size; idx4++)
	{
	  ss << idx1 << "\t";
	  ss << idx2 << "\t";
	  ss << idx3 << "\t";
	  ss << idx4 << "\t";
	  ss << TBME(idx1,idx2)(idx3,idx4) << "\n";
	}
  return ss.str();
}
