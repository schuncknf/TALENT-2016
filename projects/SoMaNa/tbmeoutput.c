#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
#include "harmon.h"
#include <math.h>

void tbmeprint (double *coeff, int nmesh)

    {
        
        int i,j, n,n1, n2, n3, n4, l1=4,l2=4,l3=4,l4=4;
        double s1,s2,s3,s4, **galagcoeff, **galagcoeffa, **bmcoeff,t1,t2,t3,t4,res;   
       	
	n = (int)coeff[0]; printf("N = %d \n\n", n);	
	
        galagcoeff = (double**) malloc(nmesh*sizeof(double*));
        for (i=0; i<nmesh; i++) { 
		galagcoeff[i] = (double*) malloc(nmesh*sizeof(double));
		}
	
	galagcoeffa = (double**) malloc((nmesh+2)*sizeof(double));
        for (i=0; i<nmesh+2; i++) { 
		galagcoeffa[i] = (double*) malloc(2*sizeof(double));
		}

	bmcoeff = (double**) malloc((nmesh+2)*sizeof(double));
        for (i=0; i<nmesh+2; i++) { 
		bmcoeff[i] = (double*) malloc((nmesh+2)*sizeof(double));
		}
	
	bmcoeff= fac_bin_mat(nmesh+4);
	
	galagcoeff= fac_galag_a0(nmesh,bmcoeff);
		
	galagcoeffa= fac_galag_a(nmesh,2,bmcoeff);
       
        FILE *TBMEOUT, *LABELSTABLE;
        TBMEOUT = fopen("tbme.out", "w");
        LABELSTABLE = fopen("states.out", "w");
       
        fprintf (TBMEOUT, "o1 o2 o3 o4 TBME \n=======================\n");
        fprintf (LABELSTABLE, "\t\t n l 2*j 2*s 2*t_z\n");
	
	/*for (i=0; i<nmesh;i++){
		for (j=0;j<nmesh;j++){
			printf("%lf ", i,j,galagcoeff[i][j]); }
			printf("# \n");}*/
	
	for (n1=0; n1<=n; n1++){//l1+=1;
       
        for (s1=+0.5; s1>-0.6; s1-=1.0){//l1+=1;
	
	for (t1=+0.5; t1>=-0.5; t1-=1.0){l1=n1*4-2*floor(s1)-floor(t1)+1; fprintf (LABELSTABLE,"Orbit number: %d %d %d %d %d %d \n", l1, n1, 0,1,(int)(2.*s1),(int)(2.*t1));}}}
       
        for (n1=0; n1<=n; n1++){//l1+=1;
       
            for(n2=0; n2<=n; n2++){//l2+=1;
           
                for (n3=0; n3<=n; n3++){//l3+=1;
       
                    for(n4=0; n4<=n; n4++){//l4+=1;
                   
                        {
                        for (s1=+0.5; s1>-0.6; s1-=1.0){
                        for (s2=+0.5; s2>-0.6; s2-=1.0){
                        for (s3=+0.5; s3>-0.6; s3-=1.0){
                        for (s4=+0.5; s4>-0.6; s4-=1.0){
					res=galag2D(nmesh,n,n1,n2,n3,n4,coeff,s1,s2,s3,s4,galagcoeff,bmcoeff,galagcoeffa);
                       	for (t1=+0.5; t1>=-0.5; t1-=1.0){l1=n1*4-2*floor(s1)-floor(t1)+1;//fprintf (LABELSTABLE,"Orbit number: %d %d %d %d %d %d \n", l1, n1, 0,1,(int)(2.*s1),(int)(2.*t1));
			for (t2=+0.5; t2>=-0.5; t2-=1.0){l2=n2*4-2*floor(s2)-floor(t2)+1;
                        for (t3=+0.5; t3>=-0.5; t3-=1.0){l3=n3*4-2*floor(s3)-floor(t3)+1;
                        for (t4=+0.5; t4>=-0.5; t4-=1.0){l4=n4*4-2*floor(s4)-floor(t4)+1;
			//printf ("%d %d %d %d \n",l1,l2,l3,l4);
                                               
                       fprintf (TBMEOUT,"%d %d %d %d %lf \n", l1, l2,l3,l4,res);
                        // printf ("%d \t\t %d \t\t \%d \t %d \t %d \t %d \t \%d \t %d \t\n", n1,(int)(2*s1),n2,(int)(2*s2),n3,(int)(2*s3),n4,(int)(2*s4));
                       
                        }}}} }}}} }}}} } 
	
        fclose(TBMEOUT);
       
        fclose(LABELSTABLE);
	
	//free(galagcoeff); free (galagcoeffa); free (bmcoeff);
    }
