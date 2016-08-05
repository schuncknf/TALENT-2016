#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>

#define hbar 197.3269788
#define pi 4*atan(1)

// Code for the eigenfunctions of a harmonic oscillator. 

// Coefficients for H.O. wave functions.

double Anl(int n,double *coeff, double **bmcoeff)

	{
		
			
		/*printf("========++++ %d %lf \n", n,pow(2.,(double)n+2));
printf("======== %lf \n", bmcoeff[0][n]);
printf("======== %lf \n", bmcoeff[1][2*n+1]);*/
		return sqrt(   pow(2.,(double)n+2)
			*bmcoeff[0][n]/
				(sqrt(3.14)*bmcoeff[1][2*n+1])	);  //return 0.1;
	
	}

// Reduced b function.

double bred(double *coeff)
	
	{
	
		return sqrt(197.3269788*197.3269788/coeff[3]/coeff[2]);

// The second hbar comes from the fact that omega is defined in our input files as hbar*omega; so w is everywhere
// actually omega*hbar and we need an extra hbar in the numerator.
	
	}

// Radial wavefunctions for the harmonic oscillator. In every function m is mass and w is circular frequency
// This function is no longer used but could be useful in the future.

double Rnl (int n, double *coeff, double **bmcoeff, double r)

	{
	
		double rnlres, br;
				
		br = bred(coeff);	
		
		rnlres = Anl(n,coeff,bmcoeff)
			/pow(br,1.5)
			*pow(r/br,coeff[1])
			*exp(-1.0*pow((r/br),2)/2.0)
			*gsl_sf_laguerre_n(n,(n+0.5), pow(r/br,2));
		
		//printf("W. FUNCTION: %lf\t %lf \n",rnlres,Anl(n,coeff,bmcoeff));		
		
		return rnlres;
	
	}
