#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "harmon.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_laguerre.h>

double** fac_bin_mat (int);

double* fac_galag_a0 (int, double** );

double* fac_galag_a (int , int a, double** );

double integrand1r (double *,double ,double);

double integrand1s (double *,double ,double);

double integrand2 (double **, int, int,int , int , int , int , double , double , double , double , double **,double**, double , double );

double galag1D (int , int , int , int , int , double **,double , double , double, double , double , double , double **, int);

double galag2D (int , int , int , int , int , double , double **, double , double , double , double , double **, double**,double **);
