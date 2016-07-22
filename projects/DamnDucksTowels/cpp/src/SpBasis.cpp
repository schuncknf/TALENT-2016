#include <cmath>
#include <gsl/gsl_sf_laguerre.h>
#include "SpBasis.h"
#include "global.h"

double sqrtFactorial(int n)
{
  if (n <= 1)
    return 1.0;
  else
    return sqrt(n) * sqrtFactorial(n-1);
}

double sqrtDoubleFactorial(int n)
{
  if (n <= 1)
    return 1.0;
  else
    return sqrt(n) * sqrtDoubleFactorial(n-2);
}

SpBasis::SpBasis(double _omega, int _nMax, int _lMax)
{
  //Defining the quantum numbers
  type = "SpBasis";
  qNames = {"n","l","m"};
  qNumSize = qNames.size();
  
  //Initializing vectors of maximum numbers
  omega = _omega;
  nMax = _nMax;
  lMax = arma::ivec(_nMax+1);
  mMax = arma::imat(_nMax+1,_lMax+1);
  
  //Defining maximum numbers and determining basis size
  basisSize = 0;
  for (int n = 0; n <= nMax; n++)
  {
    // Here to specify lMax dependency on n
    lMax(n) = _lMax;
    for (int l = 0; l <= lMax(n); l++)
    {
      mMax(n,l) = l;
      basisSize += 2*l+1;
    }
  }
  
  //Filling the quantum numbers for each state
  qNumbers = arma::imat(basisSize,qNumSize);
  int i = 0;
  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax(n); l++)
      for (int m = -mMax(n,l); m <= mMax(n,l); m++)
      {
	qNumbers(i,0) = n;
	qNumbers(i,1) = l;
	qNumbers(i,2) = m;
	i++;
      }
  
  //Calculating N normalization coefficients
  nu = NUCLEON_MASS*omega/2/HBAR;
  N = arma::vec(basisSize);
  calcN();
}

void SpBasis::calcN()
{
  //Calculating N noramalization coefficients for each basis state
  for (int i = 0; i < basisSize; i++)
  {
    int n = qNumbers(i,0);
    int l = qNumbers(i,1);
    N(i) = pow(2*nu*nu*nu/M_PI,0.25)*pow(2,0.5*n+l+1.5)*sqrtFactorial(n)*pow(nu,0.5*l)/sqrtDoubleFactorial(2*n+2*l+1);
  }
}

void SpBasis::evalRadialWaveFunction(arma::mat &wfMatrix, arma::vec r)
{
  //Calculating wave function values for each basis state and each r point provided
  wfMatrix = arma::zeros(r.n_elem,basisSize);
  for (int i = 0; i < basisSize; i++)
  {
    int n = qNumbers(i,0);
    int l = qNumbers(i,1);
    arma::vec laguerre(size(r));
    for (unsigned int j = 0; j < r.n_elem; j++)
      laguerre(j) = gsl_sf_laguerre_n(n,l+0.5,2*nu*r(j)*r(j));
    wfMatrix.col(i) = N(i) * arma::pow(r,l) % arma::exp(-nu * arma::pow(r,2)) % laguerre;
  }
}