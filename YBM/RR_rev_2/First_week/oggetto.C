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
		double v = (energy-potential_well(x))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}

	double numerov_algorithm_woods(double energy, double f0, double f_,double x){
		double v = (energy-potential_woods(x))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
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

	double potential_woods(double x){
		return -V0/(1+exp((abs(x)-6)/0.5))+V0;
	}
/*
	double get_x_woods(){
		return 
	}*/