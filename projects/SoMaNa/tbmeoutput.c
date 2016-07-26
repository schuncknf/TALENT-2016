#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
#include "harmon.h"

void tbmeprint (double m, double w, double kappar, double kappas,double kappat, double Vr,double Vs, double Vt, int n, double l)

	{
		
		int n1, n2, n3, n4;
		
		double s1,s2,s3,s4;		
		
		FILE *TBMEOUT;
		
		TBMEOUT = fopen("tbmenew.out", "w");
		
		fprintf (TBMEOUT, "n1 \t 2s1 \t n2 \t 2s2 \t n3 \t 2s3 \t n4 \t 2s4\t TBME \n");
		
		for (n1=0; n1<=n; n1++){
		
			for(n2=0; n2<=n; n2++){
			
				for (n3=0; n3<=n; n3++){
		
					for(n4=0; n4<=n; n4++){
					
						{
						
						for (s1=-0.5; s1<0.6; s1+=1.0){
						for (s2=-0.5; s2<0.6; s2+=1.0){
						for (s3=-0.5; s3<0.6; s3+=1.0){
						for (s4=-0.5; s4<0.6; s4+=1.0){
						
							fprintf (TBMEOUT,"%d \t %d \t \%d \t %d \t %d \t %d \t \%d \t %d \t %lf \n", n1,(int)(2*s1),n2,(int)(2*s2),n3,(int)(2*s3),n4,(int)(2*s4), twodgalag(5, n1,n2,n3,n4,m,w,Vr,Vs,Vt,kappar,kappas,kappat,0.0,s1,s2,s3,s4));
						
						// printf ("%d \t %d \t \%d \t %d \t %d \t %d \t \%d \t %d \t\n", n1,(int)(2*s1),n2,(int)(2*s2),n3,(int)(2*s3),n4,(int)(2*s4));
						
						}}}}} }}}}
		
		fclose(TBMEOUT);
	
	}

/* Module to run the TBME code and output the results to a table, tbme.out. The format is:

n1	n2	n3	n4	TBME.

Again, this is not the most general implementation but it's Friday and this works well enough.*/
