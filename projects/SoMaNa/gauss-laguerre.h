#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>

double fun1 (double, int, int, int, int, double, double, double);

double int1 (int, int, int, int, double, double, double, double,double);

double fun2 (double, double, int, int, int, int, double, double,double);

double int2 (int, int, int, int, double, double,double);

double tbme (int, int, int, int, double, double);

double laguerre (int, double);

int factorial (int);

int binomial (int, int);

/*

double galag (int, double*(double, int, int, int, int), int, int, int, int, double, double,double,double,double);

double twodgalag (int, double*(double, int, int, int, int,double,double), int, int, int, int, double, double,double);

*/


double galag (int, int, int, int, int, double, double,double,double,double);

double twodgalag (int, int, int, int, int, double, double,double);
