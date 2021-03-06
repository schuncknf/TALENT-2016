#include"Esercizio.h"


	double numerov_algorithm(double energy, double f0, double f_, double x){
		double v = (energy-V(x+17))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}

	double V(double x)
	{
		if (selectFunc == 1) return infSW();
		if (selectFunc == 2) return potential_well(x);
		if (selectFunc == 3) return woodsSaxon(x);
		if (selectFunc == 4) return woodsSaxon(x) + V_so(x,0,-0.5);
	}

	// Infinite square well potential
	double infSW()
	{
		return 0;
	}

	double potential_well(double x){
		double half_box = (width_box-width)/2.;
		if(x<half_box)
			return V0;
		else if(x>=half_box&&x<=6+half_box)
			return 0;
		else if(x>6+half_box)
			return V0;
	}
/*
	double woodsSaxon(double x){
		double WS;
		int         A = nProton + nNeutron;
		double     r0 = 1.27;
		double      R = r0 * pow(A,1./3);
		double coeffV = -51 + 33*((nNeutron-nProton)/A);

		if(x>width_box/2)
		{
			WS = coeffV * (1 / ( 1 + exp( ((x-(width_box/2))-R)/a_comp )));
		}
		else 
		{
			WS = coeffV * (1 / ( 1 + exp( ((-x+(width_box/2))-R)/a_comp )));
		}
		return WS;
	}
*/	
	double woodsSaxon(double s){
		return -V0/(1+exp((abs(s)-6)/0.5))+V0;
	}

	double V_so(double r, double l, double s){
		 double Vso = 0;
		  double r0 = 1.27;
		      int A = nProton + nNeutron;
		   double R = r0 * pow(A,1./3);

		   double j = l+s;
		double V_ls = j*(j+1) - l*(l+1) - 3/4;

		Vso = V_ls * r0 * r0 * (-1/r) * ((exp((r+R)/a_comp))/(a_comp*pow((exp(R/a_comp)+exp(r/a_comp)),2)));

		return Vso;
	}
