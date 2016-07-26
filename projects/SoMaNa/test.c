#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"
#include "tbmeoutput.h"
#include "read_file.h"

#define MAX_PAR 21

int main (void)

	{
		double res, *coefficients,kappar,kappat,kappas,m,w,Vr,Vt,Vs,l;
		
		int i,j,n;
		
		coefficients = (double*) malloc (MAX_PAR*sizeof(double));
		
		/* res = int1(fun1, 1,2,3,4,1.2,3.4,0.3,0.1);
		
		printf("Res1 = %f \n\n\n", res);
	
		printf("Fun2(0.1,1.0) = %f", fun2(0.1,1.0,1,2,3,4,20,30));
		
		res = int2(fun2,1,2,3,4,0.1, 0.3);
		
		printf("Res2 = %f \n\n\n", res); */
		
		coefficients=read_file();
		
		kappar = coefficients[7]; printf("\n kappa_r= %lf", kappar);
		
		Vr = coefficients[4]; printf("\n V_r= %lf", Vr);
		
		kappas = coefficients[8]; printf("\n kappa_s= %lf", kappas);
		
		Vs = coefficients[6]; printf("\n V_s= %lf", Vs);
		
		kappat = coefficients[8]; printf("\n kappa_t= %lf", kappat);
		
		Vt = coefficients[5]; printf("\n V_t= %lf", Vt);		
		
		m = coefficients[3]; printf("\n m= %lf", m);
		
		w = coefficients[2]; printf("\n omega= %lf", w);
		
		l = coefficients[1]; printf("\n l= %lf", l);
		
		n = coefficients[0]; printf("\n n= %d", n);	
		
		for (i=0; i<30; i++)
		
			{
			printf("%lf \n",Rnl(2,0,900,20,i*0.001));
			
			} 
		
		printf("\n\n STOP");		
		
		tbmeprint(m,w/197.3269788, kappar,Vr,n,l);
		
		free(coefficients);
	
	}
