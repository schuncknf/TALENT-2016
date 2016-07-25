#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"

// Module to calculate two-body matrix elements; currently it is specifically adapted to the simple V outlined in the manuals for the course. Later we will of course modify it for the more complex case, inevitably breaking everything in the process.

// The calculation proceeds in four steps, really; first the function to be integrated in the 1D integral is calculated, then it is integrated, then it is all repeated for the 2D integral (r1, r2).

// There is a lot of passing functions as arguments. Even I don't believe it's all correct.

// 1D function.

double fun1 (double x, int n1, int n2, int n3, int n4, double r1,double r2, double V)
	
	{
	
		double res;
		
		res = V*r1*r2;
		
		printf ("fun1 = %lf (%lf, %lf, %lf)", res,V,r1,r2); 
		
		return res;
	
	}

// 1D integral.

double int1 (int n1, int n2, int n3, int n4, double r1, double r2, double m, double w,double V)
	
	{
		
		double res;
		
		res=0.5*exp((-m*(pow(r1,2)+pow(r2,2))))*galag(5, n1, n2, n3, n4,m,w,r1,r2,V);
		
		printf ("int1 = %lf %lf %lf %lf %lf", res,(-1*m*(pow(r1,2)+pow(r2,2))),r1,r2,m); 		
		
		return res;
	}

// 2D function.

double fun2 (double r1, double r2, int n1, int n2, int n3, int n4, double m, double w, double V)
	
	{
	
		double res;
		
		res= Rnl(n1,0,m,w,r1)*Rnl(n2,0,m,w,r2)*Rnl(n3,0,m,w,r1)*Rnl(n4,0,m,w,r2)*exp(-m*(pow(r1,2)+pow(r2,2)))*exp(r1)*exp(r2)*int1(n1,n2,n3,n4,r1,r2,m,w,V);
		
		
		printf ("fun2 = %lf", res); 		
		
		return res;
	
	}

// 2D integral.

double int2 (int n1, int n2, int n3, int n4, double m, double w, double V)

	{
		
		double res;
		
		res=twodgalag(5,n1,n2,n3,n4,m,w,V);
		
		printf ("int2 = %lf", res); 		
		
		return res;
	
	}

// This was used for testing purposes.

/*
double tbme (double(*pot)(double, double), int n1, int n2, int n3, int n4, double m, double w)

	{
	
		double res;
		
		res = 1;
		
		return res;
	
	}*/
