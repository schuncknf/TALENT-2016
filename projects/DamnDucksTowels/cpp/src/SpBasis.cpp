#include <cmath>
#include <gsl/gsl_sf_laguerre.h>
#include "SpBasis.h"
#include "global.h"
#include "utils.h"

SpBasis::SpBasis(double _omega, int _nMax, int _lMax) :
  Basis(std::string("SpBasis"),
        std::vector<std::string>(
{"n", "l", "m", "s"
})),
omega(_omega),
nMax(_nMax),
lMax(_nMax + 1),
mMax(_nMax + 1, _lMax + 1)
{
  //Defining maximum numbers and determining basis size
  size = 0;

  for (int n = 0; n <= nMax; n++)
  {
    // Here to specify lMax dependency on n
    lMax(n) = _lMax;

    for (int l = 0; l <= lMax(n); l++)
    {
      mMax(n, l) = l;
      size += 2 * l + 1;
    }
  }

  //Considering spin
  size *= 2;
  //Filling the quantum numbers for each state
  qNumbers = arma::imat(size, qNumSize);
  int i = 0;

  for (int s = -1; s <= 1; s += 2)
    for (int n = 0; n <= nMax; n++)
      for (int l = 0; l <= lMax(n); l++)
        for (int m = -mMax(n, l); m <= mMax(n, l); m++)
        {
          qNumbers(i, 0) = n;
          qNumbers(i, 1) = l;
          qNumbers(i, 2) = m;
          qNumbers(i, 3) = s;
          i++;
        }

  //Calculating N normalization coefficients
  nu = NUCLEON_MASS * omega / 2 / HBAR;
  N = arma::zeros<arma::vec>(size);
  calcN();
}

SpBasis::~SpBasis()
{
}

int SpBasis::deltaSpin (int idx1, int idx2)
{
  if (qNumbers(idx1, 3) == qNumbers(idx2, 3))
    return 1;
  else
    return 0;
}

void SpBasis::calcN()
{
  //Calculating N noramalization coefficients for each basis state
  for (int i = 0; i < size; i++)
  {
    int n = qNumbers(i, 0);
    int l = qNumbers(i, 1);
    N(i) = pow(2 * nu * nu * nu / M_PI, 0.25) * pow(2, 0.5 * n + l + 1.5) * sqrtFactorial(n) * pow(nu, 0.5 * l) / sqrtDoubleFactorial(2 * n + 2 * l + 1);
  }
}

void SpBasis::evalRadialWaveFunction(arma::mat &wfMatrix, arma::vec &r)
{
  //Calculating wave function values for each basis state and each r point provided
  wfMatrix = arma::zeros(r.n_elem, size);

  for (int i = 0; i < size; i++)
  {
    int n = qNumbers(i, 0);
    int l = qNumbers(i, 1);
    arma::vec laguerre(r.n_elem);

    for (unsigned int j = 0; j < r.n_elem; j++)
      laguerre(j) = gsl_sf_laguerre_n(n, l + 0.5, 2 * nu * r(j) * r(j));

    wfMatrix.col(i) = N(i) * arma::pow(r, l) % arma::exp(-nu * arma::pow(r, 2)) % laguerre;
  }
}

void SpBasis::evalDerivativeRadialWaveFunction(arma::mat &wfMatrix, arma::vec &r)
{
  //Calculating first derivative of wave function for each basis state and each r point provided
  wfMatrix = arma::zeros(r.n_elem, size);

  for (int i = 0; i < size; i++)
  {
    int n = qNumbers(i, 0);
    int l = qNumbers(i, 1);
    arma::vec laguerre(r.n_elem);
    arma::vec laguerre_der(r.n_elem);

    for (unsigned int j = 0; j < r.n_elem; j++)
    {
      laguerre(j) = gsl_sf_laguerre_n(n, l + 0.5, 2 * nu * r(j) * r(j));

      if (n == 0)
        laguerre_der(j) = 0.0;
      else
        laguerre_der(j) = gsl_sf_laguerre_n(n - 1, l + 1.5, 2 * nu * r(j) * r(j));
    }

    wfMatrix.col(i) = N(i) * arma::exp(-nu * arma::pow(r, 2)) % (laguerre % (l * arma::pow(r, l - 1) - 2 * nu * arma::pow(r, l + 1)) + laguerre_der % (-4 * nu * arma::pow(r, l + 1)));
  }
}

void SpBasis::evalRadialWaveFunctionNoExp(arma::mat &wfMatrix, arma::vec &r)
{
  //Calculating wave function values for each basis state and each r point provided
  wfMatrix = arma::zeros(r.n_elem, size);

  for (int i = 0; i < size; i++)
  {
    int n = qNumbers(i, 0);
    int l = qNumbers(i, 1);
    arma::vec laguerre(r.n_elem);

    for (unsigned int j = 0; j < r.n_elem; j++)
      laguerre(j) = gsl_sf_laguerre_n(n, l + 0.5, 2 * nu * r(j) * r(j));

    wfMatrix.col(i) = N(i) * arma::pow(r, l) % laguerre;
  }
}
