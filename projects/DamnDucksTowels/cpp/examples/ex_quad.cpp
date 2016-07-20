#include <iostream>
#include <iomanip>

#include <armadillo>
#include "quadrature.h"

int main(void)
{
  arma::vec wi_p;
  arma::vec pi_p;
  arma::vec fun; 

  for (int i = 5; i < 180; i += 5)
  {
    GET_LAG_ROOTS(i,pi_p,wi_p);

    // We will integrate x \mapsto e^{-x2}(1.4x^3 + 0.3x^2 + 2.8x + 0.4)dx between 0 and +inf with gauss-laguerre:
    // calculation of the "polynomial part" of the function we want to integrate
    fun = arma::exp(pi_p) % arma::exp(-arma::square(pi_p))
      % (0.4 * arma::pow(pi_p, 3) + 0.3 * arma::pow(pi_p, 2) + 2.8 * arma::pow(pi_p, 1) + 0.4 * arma::pow(pi_p, 0));

    // This display the value of the integrale between 0 and +inf
    printf("%.3i %.15f\n", i, arma::accu(wi_p % fun));
  }
}
