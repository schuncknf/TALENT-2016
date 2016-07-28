#include"Esercizio.h"


	double numerov_algorithm_woods_proton(double energy, double f0, double f_,double r, double S, double L, double J){
		double v = (energy-potential_woods(r)-centrifug_term(r,L)-potential_spin_orbit(r,S,L,J)-potential_coulomb(r))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}


	double potential_woods(double r){
		double V = -51.+33.*(N-Z)/(N+Z); 
		double R = r0*pow(A,1./3.);
		return V/(1+exp((r-R)/0.67));
	}

	double centrifug_term(double r, double L){
		return m_factor*L*(L+1)/pow(r,2);
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
			double V = -51.+33.*(N-Z)/(N+Z); 
			double R = r0*pow(A,1./3.);
			double S_L = 1./2.*(J*(J+1)-L*(L+1)-3./4.);
			double df = -1./(pow(1+exp((r-R)/0.67),2))*(1./0.67)*exp((r-R)/0.67);
			return -0.22*V*1./r*df*S_L*pow(r0,2);
		}
		else
			return 0;
	}

	double normalise(double eigen, double n_step_width_box, double S, double L, double J){
	vector<double> wfWork;
	wfWork.push_back(0.1);
	wfWork.push_back(0.15);

	double sum = 0.;
	for(int i=1; i<n_step_width_box+1; i++) wfWork.push_back(numerov_algorithm_woods_proton(eigen, wfWork[i], wfWork[i-1], i*h_width + 2*h_width, S, L, J));
	for(int i=1; i<n_step_width_box+1; i++) sum += 4*M_PI*h_width*pow(wfWork[i],2);

	wfWork.clear();
	return sum;
}