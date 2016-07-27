#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "diag.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>
//================================ Global Constants Declareration ================================
#define N 10 //number of states in  H space
#define P 2 //number of nucleons
#define QN 4//number of quantum numbers needed for Two Body Matrix Element input
#define omega 1.0 //1*pow(10,34) // frequency of Harmonic Oscillator basis
#define maxIter 20 // maximum number of iteration for HF solver loop

//======================================== Function Declaration ==================================
double h_0(int n, int l); // function for one body hamiltonian
//double TBME(); // function from Two Body Matrix Element module
double deltaF(int a,int b);
double density(int beta, int delta, gsl_matrix_complex * D, int p);


//=========================================== Main Code ==========================================
int main()
{
  int i,j,alpha, gamma, beta, delta,n,l;
  double eval_R, evalPrev_R;
  long int row_TBME = pow(N+1,4);
  double cCrit = 1.0*pow(10.0,-5.0); //convergence criteria
  FILE * fp;
  FILE * fp0;
  fp = fopen("results.txt","wr");
  char line[128];
  char *token;
  int qn1,qn2,qn3,qn4;
  qn1 = 0; qn2 = 0; qn3 = 0; qn4 = 0;
  //============================= MATRIX & VECTOR INITIALIZATION =================================
  double TBME[21][21][21][21];//rewritten TBME into 4 dimentional matrix depending on quantum numbers, use N+1
  gsl_matrix * H = gsl_matrix_calloc (N,N); // full hamiltoninan density
  gsl_vector_complex *eval = gsl_vector_complex_alloc (N); // vector for eigenvalue
  gsl_vector_complex *evalPrev = gsl_vector_complex_alloc (N);  // vector to store previous eigenvalues
  gsl_vector * diff = gsl_vector_calloc (N); // vector for checking convergence criteria
  // initialize convergence vector elements to 1, so it won't stop loop at first run
  gsl_vector_set_all(diff, 1.0);
  //printf("diff = %f", gsl_vector_get(diff,N-1));
  gsl_matrix_complex *D = gsl_matrix_complex_alloc (N, N); // density matrix / eigenvector storeage
  gsl_matrix_complex_set_identity(D);//initialization for density matrix
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(N); //gsl diag. function workspace
  
  //============================= **END MATRIX & ARRAY INITIALIZATION** ============================
  
  //============== READ in TBME[alpha][beta][gamma][delta] matrix from .dat file ===============
  fp0 = fopen("tbme.out","r");
  i = 0; j = 0;
  //printf("row = %ld", row_TBME);
  while (fgets(line, sizeof(line), fp0))
  {
    token = strtok(line," ");
    while (token != NULL && strcmp(token,"\n"))
    {
        if (j==0){
            qn1 = strtol(token, NULL, 10);
        }
        else if (j==1){
            qn2 = strtol(token, NULL, 10);
        }
        else if (j==2){
            qn3 = strtol(token, NULL, 10);
        }
        else if (j==3){
            qn4 = strtol(token, NULL, 10);
        }
        else if (j==4){
            TBME[qn1][qn2][qn3][qn4] = strtod(token, NULL);
            //if (qn1<4 && qn2<4 && qn3<4 && qn4<4)
                //printf("TBME[%d][%d][%d][%d] = %f\n",qn1,qn2,qn3,qn4 ,TBME[qn1][qn2][qn3][qn4]);
        }
        
        token = strtok(NULL," \t ");
        j++;
        if (j == QN+1){
            j = 0;
            i++;
        }
    }
  }
  //printf("First TBME[0][3][0][0] = %f",TBME[0][3][0][0]);
  fclose(fp0);
  //============== **END READ in TBME[alpha][beta][gamma][delta] matrix from .dat file** ===============
  
  
  //==================================== HF SOLVER LOOP ==========================================
  int iter = 1;
  double spPot;
  double h0;
  int N0 = 10; // N0 is used for testing purpose, current TBME data is corrupt.
  while(iter < maxIter){
    for (alpha=0; alpha < N0; alpha++){
        for (gamma=0; gamma < N0; gamma++){
            n = alpha; //!!!need to have a n,l map with alpha gamma etc. For s-wave, l = 0, we can use n = alpha.
            l = 0;
            spPot = 0;
            h0 = deltaF(alpha, gamma) * h_0(n, l);
            for (beta=0; beta < N0; beta++){
                for (delta=0; delta < N0; delta++){
                    spPot += density(beta, delta, D, P) * TBME[alpha][beta][gamma][delta];
                    //There could be a factor of 1/2 or 1/4 etc in the spPot expression.
                    //printf("density[%d][%d] = %f\n" , beta,delta, density(beta, delta, D, P));
                    //printf("TBME[%d][%d][%d][%d] = %f\n",alpha,beta, gamma,delta,TBME[alpha][beta][gamma][delta]);
                }
            }
            gsl_matrix_set (H, alpha, gamma, (h0 + spPot));
            
        }
    }
    
    //================================== Matrix Diagonalization ==================================
    //loops for checking Hamiltoninan m.e.
    printf("Hamiltonian: \n");
    for (i = 0; i< N; i++){
        for (j=0; j<N;j++){
            double hold0;
            hold0 = gsl_matrix_get (H,i,j);
            printf("%f \t", hold0);
        }
        printf("\n");
    }

    gsl_eigen_nonsymmv(H, eval, D, w);
    gsl_eigen_nonsymmv_sort(eval, D, GSL_EIGEN_SORT_ABS_ASC);
    
    //loops for checking eigenvalues and eigenvectors.
    for (i = 0; i< N; i++){
        printf("E[%d] = %f\n  ",i, GSL_REAL(gsl_vector_complex_get(eval,i)));
        for (j=0; j<N;j++){
            gsl_complex hold;
            hold = gsl_matrix_complex_get (D,i,j);
            //printf("D[%d][%d] = %f + i* %f \n",i,j, GSL_REAL(hold), GSL_IMAG(hold));
        }
    }

   
    //=============================== **End Matrix Diagonalization** =============================
    
    //====================================== Check Convergence ===================================
    printf("Iteration: %d,   ",iter);
    for (i = 0; i < N ; i++){
        eval_R = GSL_REAL(gsl_vector_complex_get(eval,i));
        evalPrev_R = GSL_REAL(gsl_vector_complex_get(evalPrev,i));
        gsl_vector_set(diff, i, fabs(eval_R - evalPrev_R));
    }
    gsl_sort_vector(diff);
    printf("Largest diff = %f, cCrit = %f\n",gsl_vector_get(diff,N-1),cCrit);
    if ( (float)gsl_vector_get(diff,N-1) < (float)cCrit ){
        printf("Eigenvalues converge, ending calculation...\n");
        break;
    }
    //================================== ** End Check Convergence** ================================
    
    //============================================ OUTPUT ==========================================
    else {
        //printf("Iteration: %i, diff = %f *********************\n", iter, gsl_vector_get(diff,N-1));
    }
    gsl_vector_complex_memcpy (evalPrev, eval);
    //write eigenvalues to .dat file
    //gsl_vector_complex_fprintf(fp, eval, "%f");
    
    iter++;
  }
  
    //====================================== **END OUTPUT** ==========================================
    
  //==================================== **END HF SOLVER LOOP** ======================================
  gsl_eigen_nonsymmv_free(w);
  gsl_matrix_free(H); gsl_matrix_complex_free(D);
  gsl_vector_complex_free(eval); gsl_vector_complex_free(evalPrev); gsl_vector_free(diff);
  //printf("%f",eval);
  fclose(fp);
  return(0);
}
//====================================== **END Main Code** ==========================================

//=========================================== FUNCTIONS ==========================================


//Function for one body hamiltoninan (diagonal m.e. only)
double h_0(int n, int l){
    double h0 = 0;
    double hbarOmega = 10.0;//unit in MeV
    h0 = (2*n + l + 1.5) * hbarOmega
    return h0;
}


//Kronecker delta function
double deltaF(int a, int b){
    if (a == b){
        return 1.;
    }
    else{
        return 0.;
    }
}

//Density matrix
double density(int beta, int delta, gsl_matrix_complex * D, int p){
    double rho0 = 0.;
    double bb=0;
    for(int i = 0; i < P; i++){
        bb = GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get (D,beta,i), gsl_complex_conjugate(gsl_matrix_complex_get (D,delta,i))));
        //printf("D = %f + i%f\n",GSL_REAL(gsl_matrix_complex_get (D,beta,i)), GSL_IMAG(gsl_matrix_complex_get(D,beta,i)));
        //printf("D = %f + i%f\n",GSL_REAL(gsl_matrix_complex_get (D,delta,i)), GSL_IMAG(gsl_matrix_complex_get(D,delta,i)));
        rho0 = rho0 + bb;
        //printf("density = %f + i%f\n",GSL_REAL(bb), GSL_IMAG(bb));
    }
    //printf("density = %f + i%f\n",GSL_REAL(rho0), GSL_IMAG(rho0));
    return rho0;
}
//=========================================== **END FUNCTIONS** ==========================================


