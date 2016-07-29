#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "harmon.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_laguerre.h>

double** fac_bin_mat (int n)

	{
	
		double **res;
		int i,j;
		res = (double**)malloc((n+2)*sizeof(double));
		for (i=0;i<n+2;i++) {
			res[i]=(double*)malloc((n+2)*sizeof(double));
		}
		
		res[0][0] = 1.; // Initialisation of factorial array.
		res[1][0] = 1.; res[1][1] = 1.; res[1][2] = 2.;// Initialisation for the double factorials.
		
		for (i=1;i<=n;i++) {
			res[0][i]=i*res[0][i-1];
		}
		
		for (i=1;i<=n;i++) {
			res[1][i]=i*res[1][i-2];
		}
		
		for (i=2;i<n+2;i++){
		for (j=0;j<n;j++){
		
		if (i>=j+2)
			res[i][j]=res[0][i-2]/res[0][j]/res[0][i-2-j];	
		else res[i][j]=0;		
		}}
		
		return res;
		
		free(res);
	
	}

double* fac_galag_a0 (int n, double** bmcoeff)
	
	{
	
	int i;	
	
	double *coefficients, *solutions, *ws, **res;
	
	res = (double**) malloc((n+1)*sizeof(double));
		
		for (i=0; i<n+1; i++) { res[i] = (double*) malloc(2*sizeof(double));
		}
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	ws = (double*) malloc ((n)*sizeof(double)); 

	solutions = (double*) malloc ((2*(n+1))*sizeof(double));
	
	for (i=0; i<n+1; i++)
		
		{
		coefficients[i] = bmcoeff[n+2][i]*pow((-1),i)/bmcoeff[0][i]; // Coefficients for the Laguerre polynomials.
		
		printf("COEFF%d : %lf \n", i, coefficients[i]);
		}
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i++)
		
		{
			printf("SOL = %lf \n", solutions[i]);	
		}	
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = solutions[2*i]
				/(pow(n+1,2)
				*pow( gsl_sf_laguerre_n(n+1,0.,solutions[2*i]),2. )  
				); // Calculation of weights.
			
			res[0][i] = solutions[2*i]; printf ("RES 0 %lf \n", res[0][i]);
			res[1][i] = ws[i]; printf ("RES 1 %lf %lf \n", res[1][i],gsl_sf_laguerre_n(0.,n+1,(double)solutions[2*i]));		
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	return res;
	
	free(res);
	
	}

double* fac_galag_a (int n, int a, double** bmcoeff)
	
	{
	
	int i;	
	
	double *coefficients, *solutions, *ws, **res;
	
	res = (double**) malloc((n+1)*sizeof(double));
		
		for (i=0; i<n+1; i++) { res[i] = (double*) malloc(2*sizeof(double));
		}
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	ws = (double*) malloc ((n)*sizeof(double)); 

	solutions = (double*) malloc ((2*(n+1))*sizeof(double));
	
	for (i=0; i<n+1; i++)
		
		{
		coefficients[i] = bmcoeff[a+n+2][n-i]*pow((-1),i)*bmcoeff[0][n]/bmcoeff[0][i]; // Coefficients for the Laguerre polynomials.
		
		printf("COEFF%d : %lf \n", i, coefficients[i]);
		}
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i++)
		
		{
			printf("SOL = %lf \n", solutions[i]);	
		}	
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = bmcoeff[0][n+a]*solutions[2*i]/(pow(n+1,2)*pow( gsl_sf_laguerre_n(n+1,a,solutions[2*i]),2 ) * bmcoeff[0][n] ); // Calculation of weights.
			
			res[0][i] = solutions[2*i]; printf ("RES 0 %lf \n", res[0][i]);
			res[1][i] = ws[i]; printf ("RES 1 %lf \n", res[1][i]);		
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	return res;
	
	free(res);
	
	}

double integrand1r (double *coeff,double r1,double r2)
	
	{
	
		double res;
		
		// For the integral	 of the spatial form-factor next to (1-P_sigma).	
		res = 0.5*
			( coeff[4]*exp(-1*coeff[7]*(r1*r1+r2*r2)/(2*coeff[7]*(r1+0.01)*(r2+0.01)) ))/4;
		
		printf ("DEBUG: fun(%lf , %lf) = %lf \n",r1,r2,res);	
		
		return res;
	
	}

double integrand1s (double *coeff,double r1,double r2)
	
	{
	
		double res;
		
		// For the integral	 of the spatial form-factor next to (1-P_sigma).	
		res = 0.5*
			( coeff[6]*exp(-1*coeff[9]*(r1*r1+r2*r2)/(2*coeff[9]*(r1+0.01)*(r2+0.01)) ))/4;
		
		printf ("DEBUG: fun(%lf , %lf) = %lf \n",r1,r2,res);	
		
		return res;
	
	}

double integrand2 (double *coeff,int i, int j,int n1, int n2, int n3, int n4, double s1, double s2, double s3, double s4, double **galcoeff, double ** bmcoeff,double d1, double d2)
	
	{
	
		double res,res1, res2, dummy1, dummy2;	
		
		res1 = Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j])*Rnl(coeff,bmcoeff,galcoeff[0][i])*Rnl(coeff,bmcoeff,galcoeff[0][j]) * (
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
		);
		
		
		
		
		res = res1-res2;
		
		printf ("DEBUG: fun2(%lf , %lf) = %lf - %lf = %lf \n",galcoeff[0][i],galcoeff[0][j],res1,res2,res);
		
		return res;
	
	}

double galag1D (int n, int n1, int n2, int n3, int n4, double *coeff,double r1, double r2, double s1, double s2, double s3, double s4, double **galcoeff,double ** bmcoeff, int p)

	{
	
	int i;
	
	double res=0;
	
	for (i=0; i<n; i++)
		
		{
		
			if (p=0) {res+= galcoeff[1][i]*(integrand1r(coeff,r1,r2));
			}
			else {res+= galcoeff[1][i]*(integrand1s(coeff,r1,r2));
			} // Summation of weights with function values at mesh points.
								// printf("\n\n\n sum = %f \n", res);
		}
	
	return res;
	
	}


double galag2D (int n, int n1, int n2, int n3, int n4, double m, double *coeff, double s1, double s2, double s3, double s4, double **galcoeff, double**bmcoeff,double **galcoeffa)

	{
	
	int i, j;
	
	double res = 0, dummy1, dummy2;

	for (i=0; i<n; i++)
		
		for (j=0; j<n; j+=1)
		
			{
				{
			dummy1 = galag1D(n,n1,n2,n3,n4,coeff,galcoeff[0][i],galcoeff[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,0);
			dummy2 = galag1D(n,n1,n2,n3,n4,coeff,galcoeff[0][i],galcoeff[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,1);
			
			
		
			printf("DEBUG: dummy1 = %lf, dummy2 = %lf \n",dummy1,dummy2);			
			
			res+= galcoeffa[1][i]*galcoeffa[1][j]*
			integrand2(coeff,i,j,n1,n2,n3,n4,s1,s2,s3,s4,galcoeff,bmcoeff,dummy1,dummy2);
			
			// printf ("\n\n\n 2DGALAG \t sol1 = %lf sol2 = %lf wi=%lf wi=%lf res =%lf",solutions[2*i],solutions[2*j],wi,wj,res);
			}
		}
	
	return res;
	
	}

