#include <cmath>
#include <gsl/gsl_sf_laguerre.h>
#include "FullSpBasis.h"
#include "global.h"
#include "utils.h"

FullSpBasis::FullSpBasis(double _omega, int _nMax, int _lMax) :
  Basis(std::string("FullSpBasis"), std::vector<std::string>(
{
  "n", "l", "2j"
})),
omega(_omega),
      nMax(_nMax),
      lMax(_lMax)
{
  //Defining maximum numbers and determining basis size
  size = (nMax + 1) * ((lMax + 1) * 2 - 1);
  //Filling the quantum numbers for each state
  qNumbers = arma::imat(size, qNumSize);
  int i = 0;

  for (int n = 0; n <= nMax; n++)
    for (int l = 0; l <= lMax; l++)
      for (int _2j = 2 * l + 1; (_2j >= 2 * l - 1) && (_2j > 0); _2j -= 2)
      {
        qNumbers(i, 0) = n;
        qNumbers(i, 1) = l;
        qNumbers(i, 2) = _2j;
        i++;
      }

  //Calculating N normalization coefficients
  nu = HBAR * omega / HBAR2_2M;
  N = arma::zeros<arma::vec>(size);
  calcN();
}

FullSpBasis::~FullSpBasis()
{
}

void FullSpBasis::calcN()
{
  //Calculating N noramalization coefficients for each basis state
  for (int i = 0; i < size; i++)
  {
    int n = qNumbers(i, 0);
    int l = qNumbers(i, 1);
    N(i) = pow(2 * nu * nu * nu / M_PI, 0.25) * pow(2, 0.5 * n + l + 1.5) * sqrtFactorial(n) * pow(nu, 0.5 * l) / sqrtDoubleFactorial(2 * n + 2 * l + 1);
  }
}
/*
void FullSpBasis::evalRadialWaveFunction(arma::mat &wfMatrix, arma::vec &r)
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

void FullSpBasis::evalDerivativeRadialWaveFunction(arma::mat &wfMatrix, arma::vec &r)
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

void FullSpBasis::evalRadialWaveFunctionNoExp(arma::mat &wfMatrix, arma::vec &r)
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
*/
