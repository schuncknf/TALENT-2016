#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>

double laguerre (int, double);

int factorial (int);

int binomial (int, int);

double galag (int, double*(double, int, int, int, int), int, int, int, int);

double twodgalag (int, double*(double, int, int, int, int), int, int, int, int);
