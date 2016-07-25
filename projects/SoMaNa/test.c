#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"
// #include "tbme.h"

void main (void)

	{
		double res;
		
		int i,j;
		
		/* res = int1(fun1, 1,2,3,4,1.2,3.4,0.3,0.1);
		
		printf("Res1 = %f \n\n\n", res);
	
		printf("Fun2(0.1,1.0) = %f", fun2(0.1,1.0,1,2,3,4,20,30));
		
		res = int2(fun2,1,2,3,4,0.1, 0.3);
		
		printf("Res2 = %f \n\n\n", res); */
		
		for (i=0; i<30; i++)
		
			{
			printf("%lf \n",Rnl(2,0,900,20,i*0.001));
			
			} 
		
		printf("\n\n STOP");		
		
		tbmeprint(9.39565,10.0, 3);
	
	}
