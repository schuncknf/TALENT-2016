#include <stdio.h>
#include <math.h>
#include "gauss-laguerre.h"

#define h 1
#define pi 4*atan(1)

int doublefac (int n)
	
	{
	
		int k;
		
		k = n/2;
		
		return pow(2,k)*factorial(k);
	
	}

double Anl(int n, double l)

	{
	
		return sqrt((pow(2,n+l+2)*factorial(n)/(sqrt(pi)*doublefac(2*n+(int)(2*l)+1))));
	
	}

double bred(double m, double w)
	
	{
	
		return sqrt(h/m/w);
	
	}

double Rnl (int n, double l, double m, double w, double r)

	{
	
		double rnlres;
		
		rnlres = Anl(n,l)/pow(bred(m,w),1.5)*pow(r/bred(m,w),l)*exp(-pow((r/bred(m,w)),2)/2)*laguerre((int)(l+0.5), pow(r/bred(m,w),2));
		
		return rnlres;
	
	}
