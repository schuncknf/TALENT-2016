#include"Esercizio.h"


	double numerov_algorithm(double energy, double f0, double f_){
		double v = energy/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}

	double numerov_algorithm_finitewell(double energy, double f0, double f_,double x){
		double v = (energy-potential(x))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}

	double potential(double x){
		if(x<12)
			return V0;
		else if(x>=12&&x<=47)
			return 0;
		else if(x>47)
			return V0;
	}