#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_eigen.h>

#define size 4
#define eps 0.000000000000001

int findeigenvecs(complex double hmatrix[size][size]);

int notmain(void){

complex double hmatrixdiagdata[size][size] =
                  {{-2., 1., 0., 4.},
                   {1., -2., 0., 25.},
                   {24., -12., 3., 0.},
                   {0., 0., 0., 2.}};


    printf("\n");
    printf("First find the eigenvectors of the diagonal matrix\n"); 
    printf("1:\n");
    findeigenvecs(hmatrixdiagdata);
    return 0;
}				/* End of main */

int findeigenvecs(complex double hmatrix[size][size])
{
    int i, j;
    gsl_complex tempconvert;
    gsl_matrix_complex *hermmatrix = gsl_matrix_complex_alloc(size, size);
    gsl_matrix_complex *eigenvecs = gsl_matrix_complex_alloc(size, size);
    gsl_vector *eigenvals = gsl_vector_alloc(size);
    gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(size);

/* Translate the Hermitian matrix hmatrix into gsl form as hermmatrix */
    for (i = 0; i < size; i++)
	for (j = 0; j < size; j++) {
	    GSL_SET_COMPLEX(&tempconvert, creal(hmatrix[i][j]),
			    cimag(hmatrix[i][j]));
	    gsl_matrix_complex_set(hermmatrix, j, i, tempconvert);
	};

/* Here we go. */
    gsl_eigen_hermv(hermmatrix, eigenvals, eigenvecs, w);

/* Now print out the eigenvectors we have found. */
    for (i = 0; i < size; i++) {
        printf("eigenvector %d =\n",i+1);
	for (j = 0; j < size; j++) {
	    printf("%f+I*%f    ",
		   GSL_REAL(gsl_matrix_complex_get(eigenvecs, j, i)),
		   GSL_IMAG(gsl_matrix_complex_get(eigenvecs, j, i)));
	}
	printf("\n");
    }
/* Now print out the eigenvalues. */
    for (i = 0; i < size; i++){
        printf("eigenvalue %d = %f\n", i+1, gsl_vector_get(eigenvals,i));
    }

/* Free up, just to be good. */
    gsl_eigen_hermv_free(w);
    gsl_matrix_complex_free(eigenvecs);
    gsl_vector_free(eigenvals);
    gsl_matrix_complex_free(hermmatrix);

    return 0;
}
