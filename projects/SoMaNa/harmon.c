#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>

#define hbar 197.3269788
#define pi 4*atan(1)

// Code for the eigenfunctions of a harmonic oscillator. 

// Coefficients for H.O. wave functions.

double Anl(double *coeff, double **bmcoeff)

	{
		
		int aux1, aux2;
		aux1 = 2*coeff[0]+(int)(2*coeff[1])+1; aux2 = coeff[1];
	
		return sqrt(
			(   pow(2,coeff[0]+coeff[1]+2)
			*bmcoeff[0][aux2]/
				(sqrt(pi)*bmcoeff[1][aux1]) )		);
	
	}

// Reduced b function.

double bred(double *coeff)
	
	{
	
		return sqrt(197.32/coeff[3]/coeff[2]);
	
	}

// Radial wavefunctions for the harmonic oscillator. In every function m is mass and w is circular frequency.

double Rnl (double *coeff, double **bmcoeff, double r)

	{
	
		double rnlres;
		
		rnlres = Anl(coeff,bmcoeff)
			/pow(bred(coeff),1.5)
			*pow(r/bred(coeff),coeff[1])
			*exp(-1.0*pow((r/bred(coeff)),2)/2.0)
			*gsl_sf_laguerre_n(coeff[0],(coeff[1]+0.5), pow(r/bred(coeff),2));
		
		return rnlres;
	
	}
