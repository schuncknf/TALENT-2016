#include <stdlib.h>
#include <math.h>
#include "harmon.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_legendre.h>

double galag1 (int, int, int, int, int, double, double,double, double, double, double, double ,double, double, double, double, double, double, double, double, double**);

double galag2 (int, int, int, int, int , double , double ,double, double, double, double, double,double, double, double, double, double, double, double, double,double**);

double twodgalag (int, int, int, int, int, double, double, double, double, double,double, double , double , double, double, double, double, double,double**);

int delta (double, double);

double fun1 (double, double, double, double, double, double,double, double,double);

double fun1r (double, double, double, double, double, double,double, double,double);

double fun2 (double, double, double, double, double,double, double,int, int, int, int, double, double, double, double, double,double, double , double);

double** galagcoeff (int);
