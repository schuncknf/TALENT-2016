#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>

//================================ Global Constants/Matrix Declareration ================================
#define NStates 216 //number of states in TBME
#define P 8 //number of nucleons
#define maxIter 100 // maximum number of iteration for HF solver loop
#define N 108 // N should be half of NStates when we're dealing with neutrons only
#define hbarOmega 10.0//unit in MeV
#define cCrit 1.0*pow(10.0,-6.0)

double TBME[NStates+1][NStates+1][NStates+1][NStates+1] = {{{{0}}}};
int index_n[NStates+1], index_l[NStates+1], index_j[NStates+1], index_m[NStates+1], index_t[NStates+1];
int RInd[N];
//======================================== Function Declaration ==================================
double h_0(int i); // function for one body hamiltonian (diagonal elements)
double deltaF(int a,int b);//Kronecker delta
double density(int k, int l, gsl_matrix_complex * D);//calculating density matrix using eigenvector
double EHF(gsl_matrix_complex * D);//calculating Hartree-Fock energy once eigenvalues of Hamiltonian converges
//=========================================== Main Code ==========================================
int main()
{
  //============================= MATRIX & VARIABLES INITIALIZATION =================================
  
  
  gsl_matrix * H = gsl_matrix_calloc (N,N); // full hamiltoninan density
  gsl_vector_complex *eval = gsl_vector_complex_alloc (N); // vector for eigenvalue
  gsl_vector_complex *evalPrev = gsl_vector_complex_alloc (N);  // vector to store previous eigenvalues
  gsl_vector * diff = gsl_vector_calloc (N); // vector for checking convergence criteria
  gsl_vector_set_all(diff, 1.0);
  gsl_matrix_complex *D = gsl_matrix_complex_alloc (N, N); // density matrix / eigenvector storeage
  gsl_matrix_complex_set_identity(D);//initialization for density matrix
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(N); //gsl diag. function workspace
  char line[128];
  char *token, *ptr;
  int i,j,k,l,m,n,t,alpha, beta, delta, gamma, ret;
  double ret0, eHF, eval_R, evalPrev_R;
  FILE * fp0, *fp1;
  // Choosing Input file names, for _n5_l2, use NStates 216, N 108 P 8; _n5_l0, use NStates 24, N 12 P 2
  // For testing that uses _n_l, use NStates 80, N 40, P 8
  fp0 = fopen("VM-scheme_n5_l2.dat","r");
  fp1 = fopen("spM_n5_l2.dat","r");
   
  //============================= **END MATRIX & VARIABLESINITIALIZATION** ============================
  //============================= READ IN TBME MATRIX & ORBIT NUMBER VS QUANTUM NUMBER LIST ============================
  // read in index file
  m = 0;i = 0; j = 0 ; k = 0 ; l = 0 ; t = 0 ; n = 0 ;
  int count = 0;
  while (fgets(line, sizeof(line), fp1))
  {
    if (m>0){
    ret = 0;
    ret0 = 0;
    token = strtok(line,":");
    token = strtok(NULL,"");
    ret = (int) strtol(token, &ptr, 10);
    i = ret;
    ret = (int) strtol(ptr, &ptr, 10);
    j = ret;
    ret = (int) strtol(ptr, &ptr, 10);
    k = ret;
    ret = (int) strtol(ptr, &ptr, 10);
    l = ret;
    ret = (int) strtol(ptr, &ptr, 10);
    n = ret;
    ret = (int) strtol(ptr, &ptr, 10);
    t = ret;
    index_t[i] = t;
    index_n[i] = j;
    index_l[i] = k;
    index_j[i] = l;
    index_m[i] = n;
    if (t == 1){
    RInd[count] = i;
    count++;
    }
    }
    m++;
    
  }
  
  fclose(fp1);
  
  i = 0; j = 0 ; k = 0 ; l = 0 ;
  m = 0;
  while (fgets(line, sizeof(line), fp0))
  {
    if (m>2)
    {
    ret = 0;
    ret0 = 0;
    token = strtok(line,"\n");
    ret = strtol(token, &ptr, 10);
    i = ret;
    ret = strtol(ptr, &ptr, 10);
    j = ret;
    ret = strtol(ptr, &ptr, 10);
    k = ret;
    ret = strtol(ptr, &ptr, 10);
    l = ret;
    ret0 = strtod(ptr, NULL);
        //index_t is for isospin quantum number, we only record neutrons
        if (index_t[i] == 1 && index_t[j] == 1 && index_t[k] == 1 && index_t[l] == 1 && ret0 !=0.0)
        {
            TBME[i][j][k][l] = ret0;
        }
    }
    m++;
  }
  //printf("n,l number for last neutron: %d, %d\n", index_n[RInd[39]], index_l[RInd[39]]);
  fclose(fp0);
 //====================== **END READ IN TBME MATRIX & ORBIT NUMBER VS QUANTUM NUMBER LIST** ===================
 //                       *****************************************************************
 //                       *****************************************************************
 //                       *****************************************************************
 //============================================== HF SOLVER LOOP ==============================================
  int iter = 1;
  double spPot;
  double h0;
    while(iter < maxIter)
    {
        for (i=0; i < N; i++)
        {
            alpha = RInd[i];
            for (j=i; j < N; j++)
            {
                gamma = RInd[j];
                spPot = 0;
                for (k=0; k < N ; k++)
                {
                    beta = RInd[k];
                    for (l=0; l < N ; l++)
                    {
                        delta = RInd[l];
                        if (TBME[alpha][beta][gamma][delta] != 0.0)
                        {
                            spPot += density(k,l, D)  * TBME[alpha][beta][gamma][delta];
                        }
                    }
                }
                h0 = h_0(alpha) * deltaF(i,j);
                gsl_matrix_set (H, i, j, (h0 + spPot));
                gsl_matrix_set (H, j, i, gsl_matrix_get(H,i,j));
            }
        }
    //================================== Matrix Diagonalization ==================================

    
    gsl_eigen_nonsymmv(H, eval, D, w);
    gsl_eigen_nonsymmv_sort(eval, D, GSL_EIGEN_SORT_ABS_ASC);
    

    //=============================== **End Matrix Diagonalization** =============================
    
    //====================================== Check Convergence ===================================
    for (i = 0; i < N ; i++)
    {
        eval_R = GSL_REAL(gsl_vector_complex_get(eval,i));
        evalPrev_R = GSL_REAL(gsl_vector_complex_get(evalPrev,i));
        gsl_vector_set(diff, i, fabs(eval_R - evalPrev_R));
    }
    gsl_sort_vector(diff);
    printf("Iteration: %i, diff = %.14f *********************\n", iter, gsl_vector_get(diff,N-1));
    
    if ( (float)gsl_vector_get(diff,N-1) < (float)cCrit)
    {
        for (i = 0; i< P; i++)
        {
        printf("E[%d] = %.14f\t  ",i, GSL_REAL(gsl_vector_complex_get(eval,i)));
        printf("\n");
        }
        eHF = 0;
        printf("Eigenvalues converge, ending calculation...\n");
        eHF = EHF(D);
        printf("Hartree-Fock Energy = %f \n", eHF);
        break;
    }
    else if (iter+1 == maxIter)
    {
        printf("Convergence Criteria of %f not met after %d iterations, stopping solver. \n", cCrit, iter+1);
        break;
    }
    //================================== ** End Check Convergence** ================================
    
    //============================================ OUTPUT ==========================================
    
    gsl_vector_complex_memcpy (evalPrev, eval);
    //write eigenvalues to .dat file
    //gsl_vector_complex_fprintf(fp, eval, "%f");
    iter++;
}
  
//====================================== **END OUTPUT** ==========================================
//==================================== **END HF SOLVER LOOP** ====================================
 //             *****************************************************************
 //             *****************************************************************
 //             *****************************************************************
  
  gsl_eigen_nonsymmv_free(w);
  gsl_matrix_free(H); gsl_matrix_complex_free(D);
  gsl_vector_complex_free(eval); gsl_vector_complex_free(evalPrev); gsl_vector_free(diff);
  return(0);
}
//====================================== **END Main Code** ==========================================







//=========================================== FUNCTIONS ==========================================


//Function for one body hamiltoninan (diagonal m.e. only)
double h_0(int i)
{
    int n = index_n[i];
    int l = index_l[i];
    return (2*n + l + 1.5) * hbarOmega;
}


//Kronecker delta function
double deltaF(int a, int b)
{
    if (a == b){
        return 1.;
    }
    else{
        return 0.;
    }
}

//Density matrix
double density(int k, int l, gsl_matrix_complex * D)
{
    double rho0 = 0.;
    double bb=0;
    for(int i = 0; i < P; i++)
    {
        bb = GSL_REAL(gsl_complex_mul(gsl_matrix_complex_get (D,k,i), gsl_complex_conjugate(gsl_matrix_complex_get (D,l,i))));
        rho0 += bb;
    }
    return rho0;
}

//Final Hartree Fock energy calculation
double EHF(gsl_matrix_complex * D)
{
    double eKin, ePot;
    int n1,n2,n3,n4,i,j,k,l;
    n1=0;n2=0;n3=0;n4=0;
    eKin = 0; ePot = 0;
    for ( i = 0 ; i < N ; i++)
    {
        n1 = RInd[i];
        eKin += h_0(n1) * density(i,i,D);
        for ( j = 0 ; j < N ; j++)
        {
            n2 = RInd[j];
            for ( k = 0 ; k < N ; k++)
            {
                n3 = RInd[k];
                for ( l = 0 ; l < N ; l++)
                {
                    n4 = RInd[l];
                    if (TBME[n1][n2][n3][n4] != 0.0)
                    {
                        ePot += 0.5 * density(i,k,D) * density(j,l,D) * TBME[n1][n2][n3][n4];
                    }
                }
            }
        }
    }
    printf("eKin = %f, ePot = %f\n",eKin,ePot);
    return eKin + ePot;
}

//=========================================== **END FUNCTIONS** ==========================================


