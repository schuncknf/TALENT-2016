#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>


double fun1 (double, double,double, double);

double fun2 (double, double, double, double,double, double,int, int, int, int);
double laguerre (int, double);

int factorial (int);
double function (double);

int binomial (int, int);

double galag (int, int, int, int, int, double, double,double, double, double);
double twodgalag (int, int, int, int, int, double, double, double);
