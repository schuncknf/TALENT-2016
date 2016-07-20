#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"

double fun1 (double x, int n1, int n2, int n3, int n4)
	
	{
	
		return 1;
	
	}

double int1 (double(*pot)(double, double), int n1, int n2, int n3, int n4, double r1, double r2, double m, double w)
	
	{
	
		return 0.5*exp((-m*(pow(r1,2)+pow(r2,2))))*galag(5, fun1, 0, 0, 0, 0,m,w)/m;
	
	}

double fun2 (double r1, double r2, int n1, int n2, int n3, int n4, double m, double w)
	
	{
	
		return Rnl(n1,0,m,w,r1)*Rnl(n2,0,m,w,r2)*Rnl(n3,0,m,w,r1)*Rnl(n4,0,m,w,r2)*exp(-m(pow(r1,2)+pow(r2,2)));
	
	}

double int2 (double(*int1v)(double(*pot)(double, double)), int n1, int n2, int n3, int n4, double r1, double r2, double m, double w)

	{
	
		return twodgalag(5,fun2(r1,r2,n1,n2,n3,n4,m,w),0,0,0,0,m,w);
	
	}

double tbme (double(*pot)(double, double), int n1, int n2, int n3, int n4, double m, double w)

	{
	
		double res;
		
		res = 1;
		
		return res;
	
	}

void main (void)
	
	{
	
		return 0;
	
	}
