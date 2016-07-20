#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>

double laguerre (int n, double x)
	
	{
	
		/* if (n==0)
		
			return 1;
		
		else if (n==1)
			
			return (-x+1.);
		
		else
		
			return ( (2*n-1-x)*laguerre(n-1,x) - (n-1)*laguerre(n-1,x) )/n; */
	
		int i;
		
		double result = 0;		
		
		for (i=0; i<=n; i++)
			
			result += binomial(n, i)*pow((-1),i)/factorial(i)*pow(x,i);
		
		return result;
	}

int factorial (int n)
	
	{
	
		if (n==0)
			return 1;
		else if (n==1)
			return n;
		else return n*factorial(n-1);
	
	}

double function (double x)
	
	{
	
		return pow(x,2);
	
	}

int binomial (int n, int k)
	
	{
	
		return factorial(n)/factorial(n-k)/factorial(k);
	
	}

double galag (int n, double (*funcp)(double, int, int, int, int), int n1, int n2, int n3, int n4, double m, double omega)

	{
	
	int i;
	
	double *coefficients;
	
	double *solutions; 
	
	double *w;
	
	double res = 0;
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	w = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*n)*sizeof(double));
	
	for (i=0; i<=n; i++)
		
		coefficients[i] = binomial(n, i)*pow((-1),i)/factorial(i);
	
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions);
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i+=2)
	
		printf ("%f \n", solutions[i]);
	
	for (i=0; i<n; i++)
		
		{
			w[i] = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2));
						
			printf("\n\n\n %f %f \n", solutions[2*i], pow(laguerre(n+1,solutions[2*i]),2));

			printf("%f \n", w[i]);
			
			res += w[i]*(*funcp)(solutions[2*i], n1, n2, n3, n4);
		}
	
	free(w);
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res);
	
	return res;
	
	}

double twodgalag (int n, double (*funcp)(double, double, int, int, int, int), int n1, int n2, int n3, int n4, double m, double w)

	{
	
	int i, j;
	
	double *coefficients;
	
	double *solutions; 
	
	double wi, wj;
	
	double res = 0;
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*n)*sizeof(double));
	
	for (i=0; i<=n; i++)
		
		coefficients[i] = binomial(n, i)*pow((-1),i)/factorial(i);
	
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions);
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i+=2)
	
		printf ("%f \n", solutions[i]);
	
	for (i=0; i<n; i++)
		
		for (j=0; j<2*n; j+=2)
		
			{
				wi = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2));
				
				wj = solutions[2*j]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*j]),2));
		
			res += wi*wj*(*funcp)(solutions[2*i], solutions[2*j], n1, n2, n3, n4);
			}
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res);
	
	return res;
	
	}
