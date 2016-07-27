#include"Esercizio.h"

	double numerov_algorithm_woods(double energy, double f0, double f_,double r, double S, double L, double J){
		double v = (energy-potential_woods(r)-potential_spin_orbit(r,S,L,J))/m_factor-centrifug_term(r,L);
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}


	double potential_woods(double r){
		double V = -51+33*(N-Z)/A; 
		double R = r0*pow(A,1./3.);
		return V/(1+exp((r-R)/0.67));
	}

	double centrifug_term(double r, double L){
		return L*(L+1)/pow(r,2);
	}

	double potential_coulomb(double r){
		double R = r0*pow(A,1./3.);
		if(r<=R)
			return Z*pow(e,2)/(2.*R)*(3-pow(r/R,2));
		else
			return Z*pow(e,2)/r;
	}

	double potential_spin_orbit(double r, double S, double L, double J){
		if(L!=0){
			double R = r0*pow(A,1./3.);
			double S_L = J*(J+1)-L*(L-1)-3./4.;
			double df = -1./pow(1+exp((r-R)/0.67),2)*1./0.67*exp((r-R)/0.67);
			return 0.44*S_L*pow(r0,2)*1./r*df;
		}
		else
			return 0;
	}