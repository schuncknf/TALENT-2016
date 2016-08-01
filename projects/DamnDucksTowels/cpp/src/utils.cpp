#include "utils.h"
#include <cmath>

double sqrtFactorial(int n)
{
  if (n <= 1)
    return 1.0;
  else
    return sqrt(n) * sqrtFactorial(n - 1);
}

double sqrtDoubleFactorial(int n)
{
  if (n <= 1)
    return 1.0;
  else
    return sqrt(n) * sqrtDoubleFactorial(n - 2);
}
