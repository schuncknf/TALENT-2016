#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"
#include "tbmeoutput.h"
#include "read_file.h"
#include <unistd.h>

#define MAX_PAR 21
//#define NMESH 20

int main (int argc, char **argv)

	{
		double res, *coefficients;
		int i,j,n, flags, opt,nmesh;
		
    	while ((opt = getopt(argc, argv, "n:")) != -1) {
        switch (opt) {
        case 'n':
            nmesh = atof(optarg);
            break;
        default:
            fprintf(stderr, "Usage: %s [-n number of mesh points]\n", argv[0]);
            exit(EXIT_FAILURE);
        }	}	
		
		
		
		coefficients = (double*) malloc (MAX_PAR*sizeof(double));
		coefficients=read_file();
		
		tbmeprint(coefficients, nmesh);
		
		free(coefficients);
		
		return 0;
	
	}
