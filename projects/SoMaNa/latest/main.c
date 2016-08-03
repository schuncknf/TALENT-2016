#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"
#include "tbmeoutput.h"
#include "read_file.h"

#define MAX_PAR 21
#define NMESH 15

int main (void)

	{
		double res, *coefficients;
		int i,j,n;
		
		coefficients = (double*) malloc (MAX_PAR*sizeof(double));
		coefficients=read_file();
		
		tbmeprint(coefficients, NMESH);
		
		free(coefficients);
		
		return 0;
	
	}
