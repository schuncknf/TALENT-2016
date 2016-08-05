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

// Throughout the program, printf functions have been commented out; there were put into the code to facilitate debugging, allo
// -wing us to track the values calculated for every case.

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
			res[0][i]=i*res[0][i-1]; // Filling the array of factorials.
		}
		
		for (i=4;i<=n;i++) {
			res[1][i]=i*res[1][i-2]; // Filling the array of double factorials.
		}
		
		for (i=2;i<n+2;i++){
		for (j=0;j<n;j++){
		
		if (i>=j+2)
			res[i][j]=res[0][i-2]/res[0][j]/res[0][i-2-j]; // Filling the array of binomial coefficients.	
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
	
	/*for (i=0; i<2*n; i++)
		
		{
			//printf("SOL = %lf \n", solutions[i]);	
		}	*/
	
	for (i=0; i<n; i++)
		
		{
			res[0][i] = solutions[2*i];//printf ("RES 0 %lf \n", res[0][i]); // The solutions are complex numbers
// as GSL does not have an equivalent function for real polynomials. This simply takes the real part (the complex part is al
// -ways 0.
					
		}
	
	for (i=0; i<n; i++)
		
		{
			ws[i] = solutions[2*i]
				/(pow(n+1,2)
				*pow( gsl_sf_laguerre_n(n+1,0.,solutions[2*i]),2. )  
				); // Calculation of weights.
			res[1][i] = ws[i];//printf ("RES 1 %lf \n", res[1][i]);		
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
		coefficients[i] = bmcoeff[(int)(a-0.5)+n+2][n-i]*pow((-1),i)*bmcoeff[0][n]/bmcoeff[0][i]; // Coefficients for the Laguerre polynomials.
		
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
			
			res[0][i] = solutions[2*i];// printf ("RES 0 %lf \n", res[0][i]);	
		}
	for (i=0; i<n; i++)
		
		{
			ws[i] =  sqrt(4*atan(1))*bmcoeff[1][2*i+1]/pow(2,i+1)*solutions[2*i]/( pow(i+1,2)*pow( gsl_sf_laguerre_n(i+1,a,solutions[2*i]),2 ) * bmcoeff[0][i] ); // Calculation of weights.
			
			res[1][i] = ws[i]; //printf ("RES 1 %lf \n", res[1][i]);		
		}
	
	free(ws);
	
	free(solutions);
	
	free(coefficients);
	
	return res;
	
	free (res);
	}

// As of the latest version of the code these functions are no longer used but are left in the code as they could be useful with some reworking.

// BEGIN OBSOLETE

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

// END OBSOLETE

// The problem of course is that the GL quadrature is a bad choice for evaluating this integral given the boundaries. And of course it can all be done analytically.

// The integrand of the 2D integral over r1, r2. (See above for reference.)

double integrand2 (double *coeff,int i, int j,int n1, int n2, int n3, int n4, double s1, double s2, double s3, double s4, double **galcoeff, double ** bmcoeff,double d1, double d2, double ** galcoeffa)
	
	{
	
		double res,res1, res2, dummy1, dummy2, br,r11,r12,r21,r22, precoeff11, precoeff12,rs11,rs21,rs12,rs22;	
		
// Obsolete versions of the expressions:		
		
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
		
		br = bred(coeff); // Reduced length for the oscillator.

// Now the tedious part: the expressions given below were obtained by taking the Minesotta potential, multiplying it by explicit expressions for the wave functions and re-
// -arranging the expressions so the part under the integral can be reduced to f(x) x^a e-x, which can then be integrated using Gauss-Laguerre quadrature of a-th order.

// For our problem this means substituting r^2(kappa + 1/br^2) with a new variable, say x. Making the substitutions we get an expression with a=0. (?)

// Now we want to make the substitution in reverse: from the values of mesh points we need to figure out r^2/b^2 that is the argument of the Laguerre polynomials (which we
// have lazily taken from GSL.
		
		r11 = galcoeffa[0][i]/(coeff[7]*br*br+1.) ;// printf("======== %lf \n", r11);	
		r12 = galcoeffa[0][j]/(coeff[9]*br*br+1.) ;// printf("======== %lf \n", r12);
		r21 = galcoeffa[0][i]/(coeff[7]*br*br+1.) ;// printf("======== %lf \n", r21);
		r22 = galcoeffa[0][j]/(coeff[9]*br*br+1.) ;// printf("======== %lf \n", r22);
		
// These are the numerical factors which stand in front of everything; I wouldn't be terribly surprised if we missed a factor of 1/2 or 2 but fortunately or unfortunately this 
// does not affect the final results much.

		precoeff11 = d1 * 0.25*pow((1./(coeff[7]+1./(br*br))),2.)/(pow(br,6.));
		precoeff12 = d2 * 0.25*pow((1./(coeff[9]+1./(br*br))),2.)/(pow(br,6.));
		
// Now we make the reverse substitution for the argument to sinh (which we get from the initial l=0 integral over cos(theta).

		rs11 = sqrt(galcoeffa[0][i]/(coeff[7]+1./(br*br)));
		rs21 = sqrt(galcoeffa[0][j]/(coeff[7]+1./(br*br)));
		rs12 = sqrt(galcoeffa[0][i]/(coeff[9]+1./(br*br)));
		rs22 = sqrt(galcoeffa[0][j]/(coeff[9]+1./(br*br)));

// Finally, the expressions. The first is V_D, then V_E. Each is separated into terms for kappa_s and kappa_r.

		res1 = (precoeff11*sinh(2.*coeff[7]*rs11*rs21) * Anl(n1,coeff,bmcoeff)*gsl_sf_laguerre_n(n1,0.5, r11)		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,0.5, r21)
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,0.5, r11)
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,0.5, r21)
		+ precoeff12*sinh(2.*coeff[9]*rs12*rs22) *Anl(n1,coeff,bmcoeff)*gsl_sf_laguerre_n(n1,0.5, r12)		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,0.5, r22) 
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,0.5, r12)
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,0.5, r22) )*(1.-0.5*(1+s1*s2));
		
		res2 = (precoeff11*sinh(2.*coeff[7]*rs11*rs21) * Anl(n1,coeff,bmcoeff)*gsl_sf_laguerre_n(n1,0.5, r11)		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,0.5, r21)
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,0.5, r11)
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,0.5, r21)
		+ precoeff12*sinh(2.*coeff[9]*rs12*rs22) *Anl(n1,coeff,bmcoeff)*gsl_sf_laguerre_n(n1,0.5, r12)		
		* Anl(n2,coeff,bmcoeff)*gsl_sf_laguerre_n(n2,0.5, r22)
		* Anl(n3,coeff,bmcoeff)*gsl_sf_laguerre_n(n3,0.5, r22)
		* Anl(n4,coeff,bmcoeff)*gsl_sf_laguerre_n(n4,0.5, r12) )*(1.-0.5*(1+s1*s2));
			
		
		res = (res1+res2)/2; //printf("************************* %lf *** %lf *** %lf *** %lf \n",r11,r21,sinh(r11*r21), res);

// I really have no idea about the numeric prefactor here. Naively I would expect 4*pi but that makes the solution collapse entirely.
		
		//printf ("SPECIAL DEBUG!: %lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t \n",Rnl(n4,coeff,bmcoeff,galcoeff[0][i]),Rnl(n1,coeff,bmcoeff,galcoeff[0][j]),Rnl(n2,coeff,bmcoeff,galcoeff[0][i]),Rnl(n3,coeff,bmcoeff,galcoeff[0][j]),d1,pow(coeff[7],-1.5),d2,pow(coeff[7],-1.5));		
		
		//printf ("DEBUG: fun2(%lf , %lf) = %lf\n",galcoeff[0][i],galcoeff[0][j],res);
		
		return res;
	
	}

// The following function is obsolete (see above).

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
	
	//printf("++++++ %d \n", nmesh);	
	
	int i, j;
	
	double res = 0, dummy1, dummy2;

	for (i=0; i<nmesh; i++)
		
		for (j=0; j<nmesh; j++)
		
			{
				{
			//dummy1 = galag1D(nmesh,n,n1,n2,n3,n4,coeff,galcoeffa[0][i],galcoeffa[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,0);
			//dummy2 = galag1D(nmesh,n,n1,n2,n3,n4,coeff,galcoeffa[0][i],galcoeffa[0][j],s1,s2,s3,s4,galcoeff,bmcoeff,1);

// Originally, the integral over cos(theta) was evaluated numerically, but that was a bad idea and I should feel bad for thinking it. It can be evaluated analytically.
// The actual result, sinh(2 kappa r1 r2)/(kappa r1 r2), is only partly below. The other part is in integrand2().
			
			dummy1 = 1.*(+0.25*
			( coeff[4]/(coeff[7]) )); //printf("££££££££ %lf %lf \n", galcoeff[1][i],galcoeff[1][j]);
			
			dummy2 = 1.*(-0.25*
			( coeff[6]/(coeff[9]) ));
			
			//dummy1 = 2.3504*coeff[4]; dummy2 = 2.3504*coeff[6];			
			
			//printf("DEBUG: dummy1 = %lf, dummy2 = %lf \t %lf \t %lf \n",dummy1,dummy2, galcoeffa[1][i],galcoeffa[1][j]);			
			
			res+= galcoeff[1][i]*galcoeff[1][j]*integrand2(coeff,i,j,n1,n2,n3,n4,s1,s2,s3,s4,galcoeff,bmcoeff,dummy1,dummy2,galcoeff);
			
			//res+= galcoeff[1][i]*galcoeff[1][j]*sqrt(galcoeff[0][i])*sqrt(galcoeff[0][j]);			
			
			//printf ("\n\n\n 2DGALAG \t res =%lf \n",res);
			}
		}
	
	return res;
	
	}

