#include"Esercizio.h"

	double numerov_algorithm_woods_spin(double energy, double f0, double f_,double r, double S, double L){
		double v = (energy-(potential_woods(r)+potential_spin_orbit(r,S,L)))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		cout<<potential_spin_orbit(r,S,L)<<endl;
		return (a[0] * f0 - a[1] * f_) / a[2];
	}


	double potential_woods(double r){
		double V = -51+33*(N-Z)/A; 
		double R = r0*pow(A,1./3.);
		return V/(1+exp((abs(r)-R)/0.67));
	}


	double potential_spin_orbit(double r, double S, double L){
		double R = r0*pow(A,1./3.);
		double J = S+L;
		double S_L = 1./2.*(J*(J+1)-L*(L-1)-S*(S-1));
		double df = 1./pow(1+exp((abs(r)-R)/0.67),2)*1./0.67*exp((abs(r)-R)/0.67);
		return 0.44*S_L*pow(r0,2)*1/abs(r)*df;
	}