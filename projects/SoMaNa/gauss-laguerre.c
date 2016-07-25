#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>


// Module to calculate two-body matrix elements; currently it is specifically adapted to the simple V outlined in the manuals for the course. Later we will of course modify it for the more complex case, inevitably breaking everything in the process.

// The calculation proceeds in four steps, really; first the function to be integrated in the 1D integral is calculated, then it is integrated, then it is all repeated for the 2D integral (r1, r2).

// There is a lot of passing functions as arguments. Even I don't believe it's all correct.

// 1D function.

// Gauss-Laguerre integration program.

// Laguerre polynomials.

double galag (int n, int n1, int n2, int n3, int n4, double m, double w,double r1, double r2, double V);

double twodgalag (int n, int n1, int n2, int n3, int n4, double m, double w, double V);

double fun1 (double kappa, double V,double r1, double r2)

	{
	
		double res;
		
		res = V*exp(-1*kappa*(r1*r1+r2*r2))/(2*(r1+0.01)*(r2+0.01))/4;		
		
		return res;
	
	}

double fun2 (double kappa, double V, double m, double w,double r1, double r2,int n1, int n2, int n3, int n4)
	
	{
	
		double res, dummy1;
		
		dummy1 = galag(5,n1,n2,n3,n4,m,w,r1,r2,V);
		
		res = dummy1*r1*r1*r2*r2*/* Rnl(n1,0,m,w,r1)*Rnl(n2,0,m,w,r2)*Rnl(n3,0,m,w,r1)*Rnl(n4,0,m,w,r2) **/exp(r1)*exp(r2);
		
		return res;
	
	}

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
			
			result+= binomial(n, i)*pow((-1),i)/factorial(i)*pow(x,i);
		
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

/* double galag (int n, double (*funcp)(double, int, int, int, int,double,double,double), int n1, int n2, int n3, int n4, double m, double w,double r1, double r2, double V) */

double galag (int n, int n1, int n2, int n3, int n4, double m, double w,double r1, double r2, double V)

	{
	
	int i;
	
	double *coefficients;
	
	double *solutions; 
	
	double *ws;
	
	double res = 0;
	
	coefficients = (double*) malloc ((n+1)*sizeof(double)); 
	
	ws = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*n)*sizeof(double));
	
	for (i=0; i<=n; i++)
		
		{
		coefficients[i] = binomial(n, i)*pow((-1),i)/factorial(i); // Coefficients for the Laguerre polynomials.
		}
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions); // Roots of the Laguerre polynomials.
	
	gsl_poly_complex_workspace_free (dummy);
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); // Calculation of weights.
						
			printf("\n\n\n root = %f weight = %f \n", solutions[2*i], ws[i]);
			
			res+= ws[i]*fun1(1.487,200,r1,r2); // Summation of weights with function values at mesh points.
			printf("\n\n\n sum = %f \n", res);
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	return res;
	
	}

// 2D G-L quadrature. Analogous to the function above, except the summation over one index is replaced by a double-loop summation.

/* double twodgalag (int n, double (*funcp)(double, double, int, int, int, int,double,double,double), int n1, int n2, int n3, int n4, double m, double w, double V) */

double twodgalag (int n, int n1, int n2, int n3, int n4, double m, double w, double V)

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
	
	for (i=0; i<n; i++)
		
		for (j=0; j<2*n; j+=2)
		
			{
				{
				wi = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); 
	// printf("wi = %f \n", wi);
				
				wj = solutions[2*j]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*j]),2));
			 
	// printf("wj = %f \n", wj);
		
			res+= wi*wj*fun2(1.487,200.01,m,w,solutions[2*i], solutions[2*j], n1, n2, n3, n4);
			
			printf ("\n\n\n 2DGALAG \t sol1 = %lf sol2 = %lf wi=%lf wi=%lf res =%lf",solutions[2*i],solutions[2*j],wi,wj,res);
			}
		}
	
	free(solutions);
	
	free(coefficients);
	 // Freeing up memory.
	
	return res;
	
	}
