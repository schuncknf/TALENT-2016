#include <stdio.h>
#include <math.h>
#include "gauss-laguerre.h"
#include <gsl/gsl_sf.h>

#define hbar 197.3269788
#define pi 4*atan(1)

// Code for the eigenfunctions of a harmonic oscillator.

// Double (even) factorial function. In a later stage of the programme we can shunt all of these functions to a separate module.

int doublefac (int n)
	
	{
	
		int k;
		
		k = n/2;
		
		return pow(2,k)*factorial(k);
	
	}

// Coefficients for H.O. wave functions.

double Anl(int n, double l)

	{
	
		return sqrt((pow(2,n+l+2)*factorial(n)/(sqrt(pi)*doublefac(2*n+(int)(2*l)+1))));
	
	}

// Reduced b function.

double bred(double m, double w)
	
	{
	
		return sqrt(hbar/m/w);
	
	}

// Radial wavefunctions for the harmonic oscillator. In every function m is mass and w is circular frequency.

double Rnl (int n, double l, double m, double w, double r)

	{
	
		double rnlres;
		
		rnlres = Anl(n,l)/pow(bred(m,w),1.5)*pow(r/bred(m,w),l)*exp(-1.0*pow((r/bred(m,w)),2)/2.0)*gsl_sf_laguerre_n(n,(l+0.5), pow(r/bred(m,w),2));
		
		
		printf("\n\n\n n = %d l = %lf m = %lf w = %lf r = %lf \n\n\n", n,l,m,w,r);
		
		printf("\n\n\n Anl = %lf bred = %lf laguerre = %lf", Anl(n,l),bred(m,w),gsl_sf_laguerre_n(n,(l+0.5), pow(r/bred(m,w),2)));
		
		/* printf ("%lf %lf %lf %lf \n", Anl(n,l), pow(r/bred(m,w),l), exp(-0.5*pow(bred(m,w),2)), laguerre((int)(l+0.5), pow(r/bred(m,w),2)) ); */ // Used to check elements of the radial wavefunction.
		
		return rnlres;
	
	}
