// read_file
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#define P_MAX 21  //max number of constants and parameters in the input file

double* read_file(void)
{
  int i;
  double *p, tmp;
  char line[256];
  p = (double*)malloc(P_MAX * sizeof(double));

  FILE *file; 

  file = fopen("parameters.dat","r");

  for(i=0; i < P_MAX; i++)
 {
    fgets(line, sizeof(line), file);
    fscanf(file,"%lf", &tmp);
    p[i] = tmp; printf("\n \n PARAMETER %lf", tmp);
    if (p[i-1] == p[i]) break; 
  }

  fclose(file);

  return p;
}
