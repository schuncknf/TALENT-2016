#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>

// Gauss-Laguerre integration program.

// Laguerre polynomials.

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

// Factorial function.

int factorial (int n)
	
	{
	
		if (n==0)
			return 1;
		else if (n==1)
			return n;
		else return n*factorial(n-1);
	
	}

// For test purposes.

double function (double x)
	
	{
	
		return pow(x,2);
	
	}

// Binomial coefficients.

int binomial (int n, int k)
	
	{
	
		return factorial(n)/factorial(n-k)/factorial(k);
	
	}

// 1D Gauss-Laguerre quadrature. The function is sadly not general; since a lot of arguments need to be passes to the function it is specifically tailored to the problem at hand.

double galag (int n, double (*funcp)(double, int, int, int, int,double,double), int n1, int n2, int n3, int n4, double m, double omega,double r1, double r2)

	{
	
	int i;
	
	double *coefficients;
	
	double *solutions; 
	
	double *w;
	
	double res = 0;
	
	printf("Calling galag. \n");
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	w = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*n)*sizeof(double));
	
	for (i=0; i<=n; i++)
		
		coefficients[i] = binomial(n, i)*pow((-1),i)/factorial(i); // Coefficients for the Laguerre polynomials.
	
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<2*n; i+=2)
	
		{
			
			printf ("%f \n", solutions[i]);
		
		}
	
	for (i=0; i<n; i++)
		
		{
			w[i] = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); // Calculation of weights.
						
			printf("\n\n\n %f %f \n", solutions[2*i], pow(laguerre(n+1,solutions[2*i]),2));

			printf("%f \n", w[i]);
			
			res += w[i]*(*funcp)(solutions[2*i], n1, n2, n3, n4,r1,r2); // Summation of weights with function values at mesh points.
		}
	
	free(w);
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res);
	
	return res;
	
	}

// 2D G-L quadrature. Analogous to the function above, except the summation over one index is replaced by a double-loop summation.

double twodgalag (int n, double (*funcp)(double, double, int, int, int, int,double,double), int n1, int n2, int n3, int n4, double m, double w)

	{
	
	int i, j;
	
	double *coefficients;
	
	double *solutions; 
	
	double wi, wj;
	
	double res = 0;
	
	printf("Calling 2Dgalag. \n");

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
				{
				wi = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); 
	printf("wi = %f \n", wi);
				
				wj = solutions[2*j]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*j]),2));
			 
	printf("wj = %f \n", wj);
		
			res += wi*wj*funcp(solutions[2*i], solutions[2*j], n1, n2, n3, n4,m,w);
		 
	printf("res_sum = %f \n", res);
			}
		}
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res); // Freeing up memory.
	
	return res;
	
	}
