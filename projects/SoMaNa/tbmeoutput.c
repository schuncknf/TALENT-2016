#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
#include "harmon.h"
#include <math.h>

void tbmeprint (double m, double w, double kappar, double kappas,double kappat, double Vr,double Vs, double Vt, int n, double l)

    {
        
        int i,j, n1, n2, n3, n4, l1=-1,l2=-1,l3=-1,l4=-1;
       
        double s1,s2,s3,s4, **galagres1, **galagres2, **coeff;   
       
       
        coeff = (double**) malloc(n*sizeof(double));
       
        for (i=0; i<n; i++)
       
            { coeff[i] = (double*) malloc(n*sizeof(double));
       
           
            }
       
        coeff = galagcoeff(n);       
       
        galagres1 = (double**) malloc(n*sizeof(double));
       
        for (i=0; i<n; i++)
       
            { galagres1[i] = (double*) malloc(n*sizeof(double));
       
           
            }
               
        galagres2 = (double**) malloc(n*sizeof(double));
       
        for (i=0; i<n; i++)
       
            { galagres2[i] = (double*) malloc(n*sizeof(double));
       
           
            }
       
        FILE *TBMEOUT, *LABELSTABLE;
       
        TBMEOUT = fopen("tbme.out", "w");
       
        LABELSTABLE = fopen("states.out", "w");
       
        fprintf (TBMEOUT, "orbit 1 \t 2 \t\t 3 \t\t 4 \t\t TBME \n");
       
        fprintf (LABELSTABLE, "orbit \t\t n \t\t l \t\t j \t\t 2*s \n");
       
        // Calculating coefficients for Gauss-Laguerre.
       
        for (i=0; i<n; i++) {
           
            for (j=0; j<n; j++) {
           
                galagres1[i][j] = 1;
           
            }}       
       
        for (n1=0; n1<=n; n1++){//l1+=1;
       
        for (s1=-0.5; s1<0.6; s1+=1.0){l1+=1;
                       
                        fprintf (LABELSTABLE,"%d \t\t %d \t\t %d \t\t %d \t\t %d \n", l1, n1, 0,0,(int)(2*s1));
       
            for(n2=0; n2<=n; n2++){//l2+=1;
           
                for (n3=0; n3<=n; n3++){//l3+=1;
       
                    for(n4=0; n4<=n; n4++){//l4+=1;
                   
                        {
                       
                       
                       
                        for (s2=-0.5; s2<0.6; s2+=1.0){l2=n2*2+((int)(floor(s2))+1);
                        for (s3=-0.5; s3<0.6; s3+=1.0){l3=n3*2+((int)(floor(s3))+1);
                        for (s4=-0.5; s4<0.6; s4+=1.0){l4=n4*2+((int)(floor(s4))+1);
                       
                       
                       
                            fprintf (TBMEOUT,"%d \t\t %d \t\t %d \t\t %d \t\t %lf \n", l1,l2,l3,l4,twodgalag(5, n1,n2,n3,n4,m,w,Vr,Vs,Vt,kappar,kappas,kappat,0.0,s1,s2,s3,s4,coeff));
                                               
                       
                        // printf ("%d \t\t %d \t\t \%d \t %d \t %d \t %d \t \%d \t %d \t\n", n1,(int)(2*s1),n2,(int)(2*s2),n3,(int)(2*s3),n4,(int)(2*s4));
                       
                        }}}}} }}}}
       
        fclose(TBMEOUT);
       
        fclose(LABELSTABLE);
   
    }

/* Module to run the TBME code and output the results to a table, tbme.out. The format is:

n1    n2    n3    n4    TBME.

Again, this is not the most general implementation but it's Friday and this works well enough.*/
