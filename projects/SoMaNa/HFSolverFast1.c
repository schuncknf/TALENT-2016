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
#include <complex.h>
#include <sys/time.h>


//================================ Global Constants/Matrix Declareration ================================
#define NStates 216 //number of states in TBME
#define P 8 //number of nucleons
#define maxIter 100 // maximum number of iteration for HF solver loop
#define N 108 // N should be half of NStates when we're dealing with neutrons only
#define hbarOmega 10.0//unit in MeV
#define cCrit 1.0*pow(10.0,-9.0)
// Choosing Input file names, for _n5_l2, use NStates=216, N=108 P=8; _n5_l0, use NStates=24, N=12 P=2
// For testing that uses _n_l, use NStates=80, N=40, P=8
#define file0 "TBME/VM-scheme_n5_l2.dat"
#define file1 "TBME/spM_n5_l2.dat"


struct timeval start_t, end_t, start_t0, end_t0;// for obtaining time difference for performance check
double TBME[N][N][N][N] = {{{{0}}}}; //Two Body Matrix Element matrix, written from input files, initialized with 0
//1D index arrays for storing corresponding n,l,m,t quantum numbers of an orbit number ID, according to spM.dat
int index_n[NStates+1], index_l[NStates+1], index_t[NStates+1];//index_j[NStates+1], index_m[NStates+1] will not be used currently
//RInd(RealIndex) is 1D arrays for storing actual orbit number ID, since we only store neutron TBMEs, we exclude t=-1 TBMEs,
//thus we want to creat a (e.g.) 20*20*20*20 matrix instead of a 40*40*40*40 matrix, RInd is used to restore Orbit number
//so that we can obtain n,l quantum numbers together with index_n[],index_l[]. IRnd is the reverse projection.
int RInd[N], IRnd[NStates+1];
//======================================== Function Declaration ==================================
double h_0(int i); // function for one body hamiltonian (diagonal elements)
double deltaF(int a,int b);//Kronecker delta
double density(int k, int l, gsl_matrix * D);//calculating density matrix using eigenvector
double EHF(gsl_matrix * D);//calculating Hartree-Fock energy once eigenvalues of Hamiltonian converges

//=========================================== Main Function ==========================================
int main()
{
  //============================= MATRIX & VARIABLES INITIALIZATION =================================
  
  
  gsl_matrix * H = gsl_matrix_calloc (N,N); // full hamiltoninan density
  gsl_vector *eval = gsl_vector_alloc (N); // vector for eigenvalue
  gsl_vector *evalPrev = gsl_vector_alloc (N);  // vector to store previous eigenvalues
  gsl_vector * diff = gsl_vector_calloc (N); // vector for checking convergence criteria
  gsl_vector_set_all(diff, 1.0); //initialize diff to 1 to avoid hitting 'converge' first iteration
  gsl_matrix *D = gsl_matrix_alloc (N, N); // eigenvector storage
  gsl_matrix_set_identity(D);//initialization of eigenvector
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(N); //gsl diagonalization function workspace
  //strings used for writing TBME matrix from .dat files
  char *token, *ptr, line[128];
  //dummy indices
  int i,j,k,l,m,n,t,ret;
  //Value holders
  double ret0, eHF;
  FILE * fp0, *fp1;
  fp0 = fopen(file0,"r");
  fp1 = fopen(file1,"r");
   
  //************************************END MATRIX & VARIABLESINITIALIZATION************************************
  
  //======================= IMPORT TBME MATRIX & ORBIT NUMBER VS QUANTUM NUMBER LIST ======================
  //ORBIT NUMBER data IMPORT
  m = 0;i = 0; j = 0 ; k = 0 ; l = 0 ; t = 0 ; n = 0 ;
  int count = 0;
  gettimeofday(&start_t0,NULL);
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
    //we do not need index_j, index_m right now
    //index_j[i] = l;
    //index_m[i] = n;
        if (t == 1)
        {
            RInd[count] = i;
            IRnd[i] = count;
            count++;
        }
    }
    m++;
    
  }
  
  fclose(fp1);
  
  i = 0; j = 0 ; k = 0 ; l = 0 ;
  m = 0;
  //TBME data IMPORT
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
        if (index_t[i] == 1 && index_t[j] == 1 && index_t[k] == 1 && index_t[l] == 1 && ret0 !=0.0)
        {
            //The correct way to write TBME is in order of i,j,k,l
            //Here I swapped k with j so later on in the HF solver for loops, there would be lesser cache misses.
            TBME[IRnd[i]][IRnd[k]][IRnd[j]][IRnd[l]] = ret0;
        }
    }
    m++;
  }
  gettimeofday(&end_t0,NULL);
  fclose(fp0);
 //**************************END READ IN TBME MATRIX & ORBIT NUMBER VS QUANTUM NUMBER LIST************************
 //                       *****************************************************************
 //                       *****************************************************************
 //                       *****************************************************************
 //============================================== HF SOLVER LOOP ==============================================
  gettimeofday(&start_t,NULL);
  int iter = 1;
  double spPot;
  double h0;
    while(iter < maxIter)
    {
        for (i=0; i < N; i++)
        {
            // j only runs from i to N, since Hamiltoninan is symmetric
            for (j=i; j < N; j++)
            {
                spPot = 0;
                for (k=0; k < N ; k++)
                {
                    for (l=0; l < N ; l++)
                    {
                        //we don't need to run the density loop if TBME is null.
                        if (TBME[i][j][k][l] != 0.0)
                        {
                            spPot += density(k,l, D) * TBME[i][j][k][l];
                        }
                    }
                }
                h0 = h_0(i) * deltaF(i,j);
                gsl_matrix_set (H, i, j, (h0 + spPot));
                //exploiting the symmetricity of the Hamiltonian, true for Real Hamiltonian
                gsl_matrix_set (H, j, i, gsl_matrix_get(H,i,j));
            }
        }
    //================================== Matrix Diagonalization ==================================
    gsl_eigen_symmv(H, eval, D, w);
    gsl_eigen_symmv_sort(eval, D, GSL_EIGEN_SORT_ABS_ASC);
    //************************************End Matrix Diagonalization************************************
    
    
    //====================================== Check Convergence ===================================
    
    //Calculate difference between current eigenvalue and previous eigenvalue, store into diff
    for (i = 0; i < N ; i++)
    {
        gsl_vector_set(diff, i, fabs(gsl_vector_get(eval,i) - gsl_vector_get(evalPrev,i)));
    }
    gsl_sort_vector(diff);
    printf("Iteration: %i, diff = %.12f *********************\n", iter, gsl_vector_get(diff,N-1));
    
    //Compare largest element in diff with cCrit.
    if ( (float)gsl_vector_get(diff,N-1) < (float)cCrit)
    {
        
        for (i = 0; i< P; i++)
        {
            printf("E[%d] = %.12f MeV\t  ",i, gsl_vector_get(eval,i));
            printf("\n");
        }
        eHF = 0;
        printf("Eigenvalues converge, ending calculation......\n");
        eHF = EHF(D);
        printf("Hartree-Fock Energy = %f MeV\n", eHF);
        printf("Convergence Criteria is: 10^%d\n", (int)log10(cCrit));
        gettimeofday(&end_t,NULL);
        break;
    }
    //End loop if iteration goes over pre-set limit.
    else if (iter+1 == maxIter)
    {
        printf("Convergence Criteria of %f not met after %d iterations, stopping solver...... \n", cCrit, iter+1);
        break;
    }
    //*************************************** End Check Convergence***************************************
    
    //Store current eigenvalue so as to compare with next set of eigenvalues.
    gsl_vector_memcpy (evalPrev, eval);
    iter++;
}
  //Recording time of HF solver loop & data import
  printf("Data import time (mainly TBME) = %f seconds\n", (end_t0.tv_sec - start_t0.tv_sec)*1.0 + (end_t0.tv_usec - start_t0.tv_usec)/1000000.0);
  printf("HFSolver Time (excluding data import time) = %f seconds\n", (end_t.tv_sec - start_t.tv_sec)*1.0 + (end_t.tv_usec - start_t.tv_usec)/1000000.0 );
  
  
 //*********************************************END HF SOLVER LOOP************************************************
 //                     *****************************************************************
 //                     *****************************************************************
 //                     *****************************************************************
  
 //free variables
  gsl_eigen_symmv_free(w);
  gsl_matrix_free(H); gsl_matrix_free(D);
  gsl_vector_free(eval); gsl_vector_free(evalPrev); gsl_vector_free(diff);
  return(0);
}
//*******************************************END Main Code*********************************************







//=========================================== FUNCTIONS ==========================================


//Function for one body hamiltoninan (diagonal m.e. only)
double h_0(int i)
{
    int n = index_n[RInd[i]];
    int l = index_l[RInd[i]];
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

//Density matrix calculated using eigenvectors
double density(int k, int l, gsl_matrix * D)
{
    double rho0 = 0.;
    double bb=0;
    for(int i = 0; i < P; i++)
    {
        bb = gsl_matrix_get(D,k,i) * gsl_matrix_get(D,l,i);
        rho0 += bb;
    }
    return rho0;
}

//Final Hartree Fock energy calculation
double EHF(gsl_matrix * D)
{
    double eKin, ePot;
    int i,j,k,l;
    eKin = 0; ePot = 0;
    for ( i = 0 ; i < N ; i++)
    {
        eKin += h_0(i) * density(i,i,D);
        for ( k = 0 ; k < N ; k++)
        {
            for ( j = 0 ; j < N ; j++)
            {
                for ( l = 0 ; l < N ; l++)
                {
                    if (TBME[i][k][j][l] != 0.0)
                    {
                        ePot += 0.5 * density(i,k,D) * density(j,l,D) * TBME[i][k][j][l];
                    }
                }
            }
        }
    }
    //printf("eKin = %f, ePot = %f\n",eKin,ePot);
    return eKin + ePot;
}

//*****************************************************END FUNCTIONS*****************************************************
//***************************************************** END PROGRAM *****************************************************
