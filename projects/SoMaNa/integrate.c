#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>

#define PI 3.14159265
#define QUAD(a) ( (a) * (a) )
#define CUBE(a) ( (a) * (a) * (a) )
#define QUAR(a) ( (a) * (a) * (a) * (a) )

double function(int, int, double, int), get_integral(double, double, int, int, int);
int factorial(int);

double gsl_sf_legendre_sphPlm(const int l, const int m, const double angle);

int main(void)
{
	double res, term1, term2, t1, t2, t3, t4, term;
	int l, m, l1, m1, l_max;
	
	printf("Insert l_max\n");
	scanf( "%d" , &l_max);
	printf("l_max = %i \n", l_max);

        
	res = 0.;
	for (l = 0; l <= l_max; l++){
		for (m = -l; m <= l; m++){	
			l1 = l;
			for (m1 = -l1; m1 <= l1; m1++){
				printf(" l=%d , m=%d \n", l, m);
				printf(" l1=%d , m1=%d \n", l1, m1);
				// phi part
				t1 = get_integral(0.,2.*PI, (double)l, (double)m, 1);
				t2 = get_integral(0.,2.*PI, (double)l, (double)m, 2);
				t3 = get_integral(0.,2.*PI, (double)l1, (double)m1, 1);
				t4 = get_integral(0.,2.*PI, (double)l1, (double)m1, 2);
				term1 = t1*t3 + t2*t4;
				//term1 = get_integral(0.,2.*PI, (double)l, (double)m, 4);
				//printf("term1 = %lf \n", term1);
				// integrals over theta
				term2 =  get_integral(-1., 1., l, m, 3)*get_integral(-1., 1., l1, m1, 3);
				//printf("term2 = %lf \n", term2);
				term = term1*term2;
				if (term != 0.0){
					printf(" \n term = %lf \n", term);
				}
				res = res + term;
			}
		}
	}

	printf("\n res = %lf \n", res);
	return 0;
}

double function(int l, int m, double angle, int fun){
  double func;

 switch (fun)
    {
    case 1: { /************ phi part ***********/
	m=(double)m;
	func = cos(m*angle);
      	break;}
      
    case 2: { /***********  phi part ***********/
    	m=(double)m;
	func = sin(fabs(m)*angle);
      	break;}

    case 3: { /*********  theta part ***********/
	func = gsl_sf_legendre_sphPlm(l, fabs(m), angle);
	//printf(" func = %lf \n", func);
      	break;}

    default:{
      if ((fun < 1) || (fun > 3)) {
	printf("\n error! ");	
	break;
      }
    }
    }
  return func;
}

double get_integral(double a, double b, int l, int m, int fun){
  double get_int,t,t_min,t_max,dt,tmp,sum,e,de;
  int i, n;

  n=10000;

  if (1 == 1) {
    /*low precision*/
    de = (b - a)/(double)n;
    sum = 0.5*(function(l,m,a,fun)+function(l,m,b,fun));
    for (i=1; i <= (n-1); i++)
      {
      e = a + de*(double)i;
      sum = sum + function(l,m,e,fun);
      }
    get_int = sum*de;
    //printf(" \n get_int = %14.10e", get_int);
  }
  else {
    //higher precision
    t_min = 0.;
    t_max = sqrt(b - a);
    dt = (t_max - t_min)/(double)n;
    tmp = t_max*function(l,m,b,fun);
    sum = 0.5*tmp;
    for (i=1; i <= (n-1); i++)
      {
	t = t_min + dt*(double)i;
	e = a + QUAD(t);
	sum = sum +t*function(l,m,e,fun);
      }
    get_int = sum*dt;
  }

  return get_int;
}

int factorial (int n)
	{
		if (n==0)
			return 1;
		else if (n==1)
			return n;
		else return n*factorial(n-1);
	}
