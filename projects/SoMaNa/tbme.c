#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"

// Module to calculate two-body matrix elements; currently it is specifically adapted to the simple V outlined in the manuals for the course. Later we will of course modify it for the more complex case, inevitably breaking everything in the process.

// The calculation proceeds in four steps, really; first the function to be integrated in the 1D integral is calculated, then it is integrated, then it is all repeated for the 2D integral (r1, r2).

// There is a lot of passing functions as arguments. Even I don't believe it's all correct.

// 1D function.

double fun1 (double x, int n1, int n2, int n3, int n4, double r1,double r2)
	
	{
	
		return r1*r2;
	
	}

// 1D integral.

double int1 (double(*pot)(double, int, int, int, int,double,double), int n1, int n2, int n3, int n4, double r1, double r2, double m, double w)
	
	{
	
		return 0.5*exp((-m*(pow(r1,2)+pow(r2,2))))*galag(5, fun1, 0, 0, 0, 0,m,w,r1,r2)/m;
	
	}

// 2D function.

double fun2 (double r1, double r2, int n1, int n2, int n3, int n4, double m, double w)
	
	{
	
		return Rnl(n1,0,m,w,r1)*Rnl(n2,0,m,w,r2)*Rnl(n3,0,m,w,r1)*Rnl(n4,0,m,w,r2)*exp(-m*(pow(r1,2)+pow(r2,2)))*exp(r1)*exp(r2);
	
	}

// 2D integral.

double int2 (double(*int1v)(double(*pot)(double, int, int, int, int),int, int, int, int, double, double, double, double), int n1, int n2, int n3, int n4, double m, double w)

	{
	
		return twodgalag(5,fun2,n1,n2,n3,n4,m,w);
	
	}

// This was used for testing purposes.

/*
double tbme (double(*pot)(double, double), int n1, int n2, int n3, int n4, double m, double w)

	{
	
		double res;
		
		res = 1;
		
		return res;
	
	}*/
