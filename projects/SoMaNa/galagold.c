#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>


// Module to calculate two-body matrix elements; currently it is specifically adapted to the simple V outlined in the manuals for the course. Later we will of course modify it for the more complex case, inevitably breaking everything in the process.

// The calculation proceeds in four steps, really; first the function to be integrated in the 1D integral is calculated, then it is integrated, then it is all repeated for the 2D integral (r1, r2).

// There is a lot of passing functions as arguments. Even I don't believe it's all correct.

// 1D function.


double galag (int, int, int, int, int, double, double,double,double,double);

double twodgalag (int, int, int, int, int, double, double,double);

double fun1 (double x, int n1, int n2, int n3, int n4, double r1,double r2, double V)
	
	{
	
		double res;
		
		res = -2*r1*r2;
		
		printf ("fun1 = %lf (%lf, %lf, %lf) \n", res,V,r1,r2); 
		
		return res;
	
	}

// 1D integral.

double int1 (int n1, int n2, int n3, int n4, double r1, double r2, double m, double w,double V)
	
	{
		
		double res;
		
		res=V*0.5*14.32305/(2*m*(r1+0.001)*(r2+0.001));
		
		printf ("int1 = %lf %lf %lf %lf \n", res, V, (2*1.486*(r1+0.001)*(r2+0.001))); 		
		
		return res;
	}

// 2D function.

double fun2 (double r1, double r2, int n1, int n2, int n3, int n4, double m, double w, double V)
	
	{
	
		double res;
		
		res= Rnl(n1,0,m,w,r1)*Rnl(n2,0,m,w,r2)*Rnl(n3,0,m,w,r1)*Rnl(n4,0,m,w,r2)*int1(n1,n2,n3,n4,r1,r2,m,w,V);
		
		printf ("fun2 = %lf \t Rnl = %lf \n", res,Rnl(n1,0,m,w,r1)); 		
		
		return res;
	
	}

// 2D integral.

double int2 (int n1, int n2, int n3, int n4, double m, double w, double V)

	{
		
		double res;
		
		res=twodgalag(9,n1,n2,n3,n4,m,w,V);
		
		printf ("int2 = %lf \n", res); 		
		
		return res;
	
	}

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

/* double galag (int n, double (*funcp)(double, int, int, int, int,double,double,double), int n1, int n2, int n3, int n4, double m, double w,double r1, double r2, double V) */

double galag (int n, int n1, int n2, int n3, int n4, double m, double w,double r1, double r2, double V)

	{
	
	int i;
	
	double *coefficients;
	
	double *solutions; 
	
	double *ws;
	
	double res = 0;
	
	printf("Calling galag. \n");
	
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
			printf("\n solutions = %f", solutions[i]);
		}
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); // Calculation of weights.
						
			//printf("\n\n\n root = %f weight = %f \n", solutions[2*i], ws[i]);
			
			res += ws[i]*1; // Summation of weights with function values at mesh points.
			//printf("\n\n\n sum = %f \n", res);
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res);
	
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
	
	printf("Calling 2Dgalag. \n");

	coefficients = (double*) malloc ((n+1)*sizeof(double)); 

	solutions = (double*) malloc ((2*n)*sizeof(double));
	
	for (i=0; i<=n; i++)
		
		coefficients[i] = binomial(n, i)*pow((-1),i)/factorial(i);
	
	gsl_poly_complex_workspace *dummy = gsl_poly_complex_workspace_alloc (n+1);
	
	gsl_poly_complex_solve (coefficients, n+1, dummy, solutions);
	
	gsl_poly_complex_workspace_free (dummy);
	
		for (i=0; i<n; i++)
		
		{
			printf("\n solutions = %f", solutions[i]);
		}
	
	for (i=0; i<n; i++)
		
		for (j=0; j<2*n; j+=2)
		
			{
				{
				wi = solutions[2*i]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*i]),2)); 
	// printf("wi = %f \n", wi);
				
				wj = solutions[2*j]/(pow(n+1, 2)*pow(laguerre(n+1,solutions[2*j]),2));
			 
	// printf("wj = %f \n", wj);
		
			res += wi*wj*fun2(solutions[2*i], solutions[2*j], n1, n2, n3, n4,m,w,V);
	
	
	printf ("root1 = %lf weight1 = %lf root2 = %lf weight2 = %lf\n", solutions[2*i],wi,solutions[2*j],wj);
	
		 
	printf("res_sum = %f \n", res);
			}
		}
	
	free(solutions);
	
	free(coefficients);
	
	printf("%f \n", res); // Freeing up memory.
	
	return res;
	
	}
