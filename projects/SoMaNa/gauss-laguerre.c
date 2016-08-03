#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "harmon.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_laguerre.h>

//
//
//======================GAUSS-LAGUERRE.C======================================================================================
//
// The module contains functions for evaluating integrals that are necessary to calculate two-body matrix elements, using Gauss
// -Laguerre quadrature.

// The following function forms a matrix containing factorials (first row), double factorials (second row) and binomial coeffic
// -ients (remainder). We use this instead of calling recursively-defined functions in order to gain on execution speed.

double** fac_bin_mat (int n)

	{
	
		double **res;
		int i,j;
		res = (double**)malloc((n+2)*sizeof(double));
		for (i=0;i<n+2;i++) {
			res[i]=(double*)malloc((n+2)*sizeof(double));
		}
		
		res[0][0] = 1.; // Initialisation of factorial array.
		res[1][0] = 1.; res[1][1] = 1.; res[1][2] = 2.; res[1][3] = 3.;// Initialisation for the double factorials.
		
		for (i=1;i<=n;i++) {
			res[0][i]=i*res[0][i-1];
		}
		
		for (i=4;i<=n;i++) {
			res[1][i]=i*res[1][i-2];
		}
		
		for (i=2;i<n+2;i++){
		for (j=0;j<n;j++){
		
		if (i>=j+2)
			res[i][j]=res[0][i-2]/res[0][j]/res[0][i-2-j];	
		else res[i][j]=0;		
		}}
		
		return res;
	
	}

// The function calculates mesh points and weights for a=0 Gauss-Laguerre quadrature.

double** fac_galag_a0 (int n, double** bmcoeff)
	
	{
	
	int i,j;	
	
	double *coefficients, *solutions, *ws;
	
	double **res;
	
	res = (double**) malloc((n)*sizeof(double*));
		
		for (i=0; i<n; i++) { res[i] = (double*) malloc(2*sizeof(double));
		}
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	ws = (double*) malloc ((n+1)*sizeof(double)); // Weights are stored here.

	solutions = (double*) malloc ((2*(n+1))*sizeof(double));
	
	for (i=0; i<n+1; i++)
		
		{
		coefficients[i] = bmcoeff[n+2][i]*pow((-1),i)/bmcoeff[0][i]; // Coefficients for the Laguerre polynomials.
		
		//printf("COEFF%d : %lf \n", i, coefficients[i]);
		}
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i++)
		
		{
			//printf("SOL = %lf \n", solutions[i]);	
		}	
	
	for (i=0; i<n; i++)
		
		{
			res[0][i] = solutions[2*i];printf ("RES 0 %lf \n", res[0][i]);
					
		}
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = solutions[2*i]
				/(pow(n+1,2)
				*pow( gsl_sf_laguerre_n(n+1,0.,solutions[2*i]),2. )  
				); // Calculation of weights.
			res[1][i] = ws[i];printf ("RES 1 %lf \n", res[1][i]);		
		}
	
	/*for (i=0; i<n;i++){
		for (j=0;j<n;j++){
			printf("%lf  ", i,j,res[i][j]); }
			printf("$ \n");}*/	
	
	return res;
	
	}

// Same as the above, but the mesh points and weights are calculated for the generalised Gauss-Laguerre quadrature.

double** fac_galag_a (int n, double a, double** bmcoeff)
	
	{
	
	int i;	
	
	double *coefficients, *solutions, *ws;
	
	static double **res;
	
	res = (double**) malloc((n+1)*sizeof(double));
		
		for (i=0; i<n+1; i++) { res[i] = (double*) malloc((n+1)*sizeof(double));
		}
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	ws = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*(n+1))*sizeof(double));
	
	for (i=0; i<n+1; i++)
		
		{
		coefficients[i] = sqrt(4*atan(1))*bmcoeff[(int)(a-0.5)+n+2][n-i]*pow((-1),i)*bmcoeff[0][n]/bmcoeff[0][i]; // Coefficients for the Laguerre polynomials.
		
		//printf("COEFF%d : %lf \n", i, coefficients[i]);
		}
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i++)
		
		{
			//printf("SOL = %lf \n", solutions[i]);	
		}	
	
	for (i=0; i<n; i++)
		
		{
			
			res[0][i] = solutions[2*i]; //printf ("RES 0 %lf \n", res[0][i]);	
		}
	for (i=0; i<n; i++)
		
		{
			ws[i] = bmcoeff[0][n+(int)(a-0.5)]*solutions[2*i]/(pow(n+1,2)*pow( gsl_sf_laguerre_n(n+1,a,solutions[2*i]),2 ) * bmcoeff[0][n] ); // Calculation of weights.
			
			res[1][i] = ws[i]; //printf ("RES 1 %lf \n", res[1][i]);		
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	printf ("AAAA! \n");
	return res;
	
	free (res);
	}

// The integrands of the first integral (over r) (see HF_truncated.pdf). First is the function for the r part, then the s part.

double integrand1r (double *coeff,double r1,double r2)
	
	{
	
		double res;
		
		// For the integral	 of the spatial form-factor next to (1-P_sigma).	
		res = 0.25*
			( coeff[4]/(coeff[7])); 
		/*res = 0.5*
			( coeff[4]*exp(-1*coeff[7]/(2*coeff[7]*(r1+0.01)*(r2+0.01)) ))/4;*/
		//printf ("DEBUG: fun(%lf , %lf) = %lf \t %lf \t %lf \n",r1,r2,res,exp(-1*coeff[7]*(r1*r1+r2*r2)),(2*coeff[7]*(r1+0.1)*(r2+0.1)));	
		
		return res;
	 
	}

double integrand1s (double *coeff,double r1,double r2)
	
	{
	
		double res;
		
		// For the integral	 of the spatial form-factor next to (1-P_sigma).	
		res = -0.25*
			( coeff[6]/(coeff[9]) );
		/*res = 0.5*
			( coeff[6]*exp(-1*coeff[9]/(2*coeff[9]*(r1+0.01)*(r2+0.01)) ))/4;*/
		//printf ("DEBUG: fun(%lf , %lf) = %lf \t %lf \t %lf  \t %lf\n",r1,r2,res,exp(-1*coeff[7]*(r1*r1+r2*r2)),(2*coeff[7]*(r1+0.1)*(r2+0.1)), coeff[7]);	

		
		return res;
	
	}

// The integrand of the 2D integral over r1, r2. (See above for reference.)

double integrand2 (double *coeff,int i, int j,int n1, int n2, int n3, int n4, double s1, double s2, double s3, double s4, double **galcoeff, double ** bmcoeff,double d1, double d2, double ** galcoeffa)
	
	{
	
		double res,res1, res2, dummy1, dummy2, br;	
		
		/*res1 = Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j])*Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j]) * (
		d1*galcoeff[0][i]*galcoeff[0][j]/coeff[7]
		* exp(-1*coeff[7]* (pow(galcoeff[0][i],2) + pow(galcoeff[0][j],2) - galcoeff[0][i] - galcoeff[0][j] )  )
		+ d2*galcoeff[0][i]*galcoeff[0][j]/coeff[7]
		* exp(-1*coeff[7]*(pow(galcoeff[0][i],2) + pow(galcoeff[0][j],2) - galcoeff[0][i] - galcoeff[0][j] )  )
		);
		
		res1 = Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j])*Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j]) * (
		d1*galcoeff[0][i]*galcoeff[0][j]/coeff[7]
		* exp(-1*coeff[7]* (pow(galcoeff[0][i],2) + pow(galcoeff[0][j],2) - galcoeff[0][i] - galcoeff[0][j] )  )
		+ d2*galcoeff[0][i]*galcoeff[0][j]/coeff[7]
		* exp(-1*coeff[7]*(pow(galcoeff[0][i],2) + pow(galcoeff[0][j],2) - galcoeff[0][i] - galcoeff[0][j] )  )
		);*/
		
		/*res1 = Rnl(n1,coeff,bmcoeff,galcoeff[0][i])*Rnl(n2,coeff,bmcoeff,galcoeff[0][j])*Rnl(n3,coeff,bmcoeff,galcoeff[0][i])*Rnl(n4,coeff,bmcoeff,galcoeff[0][j]) * (
		d1*(pow(2*coeff[7],-3.))+ d2*(pow(2*coeff[9],-3.)))*(1.-0.5*(1+s3*s4));
		
		res2 = Rnl(n1,coeff,bmcoeff,galcoeff[0][i])*Rnl(n2,coeff,bmcoeff,galcoeff[0][j])*Rnl(n4,coeff,bmcoeff,galcoeff[0][j])*Rnl(n3,coeff,bmcoeff,galcoeff[0][i]) * (
		d1*(pow(2*coeff[7],-3.))+d2*(pow(2*coeff[9],-3.)))*(1.-0.5*(1+s3*s4));*/
		
		br = bred(coeff);		
		
		res1 = (pow(br,-6.))*(d1 * pow((1/(coeff[7]+1/(2*br))),3) * pow(-2*coeff[7],2) 
		+ d2 * pow((1/(coeff[9]+1/(2*br))),3) * pow(-2*coeff[9],2) )*(1.-0.5*(1+s3*s4)) * Anl(n1,coeff,bmcoeff)*gsl_sf_laguerre_n(n1,(0+0.5), pow(galcoeffa[0][i]/br,2))		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,(0+0.5), pow(galcoeffa[0][j]/br,2))
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,(0+0.5), pow(galcoeffa[0][i]/br,2))
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,(0+0.5), pow(galcoeffa[0][i]/br,2));
		
		res2 = (pow(br,-6.))*(d1 * pow((1/(coeff[7]+1/(2*br))),3) * pow(-2*coeff[7],2) +
		+ d2 * pow((1/(coeff[9]+1/(2*br))),3) * pow(-2*coeff[9],2) * Anl(n1,coeff,bmcoeff) )*(1.-0.5*(1+s3*s4))*gsl_sf_laguerre_n(n1,(0+0.5), pow(galcoeffa[0][i]/br,2))		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,(0+0.5), pow(galcoeffa[0][j]/br,2))
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,(0+0.5), pow(galcoeffa[0][i]/br,2))
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,(0+0.5), pow(galcoeffa[0][j]/br,2));
			
		
		res = 1*(res1+res2);
		
		//printf ("SPECIAL DEBUG!: %lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t \n",Rnl(n4,coeff,bmcoeff,galcoeff[0][i]),Rnl(n1,coeff,bmcoeff,galcoeff[0][j]),Rnl(n2,coeff,bmcoeff,galcoeff[0][i]),Rnl(n3,coeff,bmcoeff,galcoeff[0][j]),d1,pow(coeff[7],-1.5),d2,pow(coeff[7],-1.5));		
		
		//printf ("DEBUG: fun2(%lf , %lf) = %lf\n",galcoeff[0][i],galcoeff[0][j],res);
		
		return res;
	
	}

// The 1D Gauss-Laguere integral. This is reduced to a simple sum of the integrand multiplied by weights, at the mesh points.

double galag1D (int nmesh,int n, int n1, int n2, int n3, int n4, double *coeff,double r1, double r2, double s1, double s2, double s3, double s4, double **galcoeff,double ** bmcoeff, int p)

	{
	
	int i;
	
	double res=0;
	
	for (i=0; i<nmesh; i++)
		
		{
			
			
			if (p=0) {res+= galcoeff[1][i]*(integrand1r(coeff,r1,r2));
			}
			else {res+= galcoeff[1][i]*(integrand1s(coeff,r1,r2));
			} // Summation of weights with function values at mesh points.
								//printf("\n\n\n %lf sum = %f \n", galcoeff[1][i],res);
		}
	
	return res;
	
	}


double galag2D (int nmesh,int n, int n1, int n2, int n3, int n4, double *coeff, double s1, double s2, double s3, double s4, double **galcoeff, double**bmcoeff,double **galcoeffa)

	{
	
	int i, j;
	
	double res = 0, dummy1, dummy2;

	for (i=0; i<nmesh; i++)
		
		for (j=0; j<nmesh; j+=1)
		
			{
				{
			//dummy1 = galag1D(nmesh,n,n1,n2,n3,n4,coeff,galcoeffa[0][i],galcoeffa[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,0);
			//dummy2 = galag1D(nmesh,n,n1,n2,n3,n4,coeff,galcoeffa[0][i],galcoeffa[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,1);
			
			dummy1 = 2.3504*(0.25*
			( coeff[4]/(coeff[7]) ));
			
			dummy2 = 2.3504*(-0.25*
			( coeff[6]/(coeff[9]) ));
		
			//printf("DEBUG: dummy1 = %lf, dummy2 = %lf \t %lf \t %lf \n",dummy1,dummy2, galcoeffa[1][i],galcoeffa[1][j]);			
			
			res+= galcoeff[1][i]*galcoeff[1][j]*
			integrand2(coeff,i,j,n1,n2,n3,n4,s1,s2,s3,s4,galcoeff,bmcoeff,dummy1,dummy2,galcoeffa );
			
			//printf ("\n\n\n 2DGALAG \t res =%lf \n",res);
			}
		}
	
	return res;
	
	}

