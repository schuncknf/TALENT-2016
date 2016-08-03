#include"Esercizio.h"


/////////////////////HARTREE FOCK///////////////////////////

	double numerov_algorithm_HF(double energy, double f0, double f_,double r, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1){
		double v=0., a[3];
			a[0] = 2. * (1. - 5./12. * v_HF(energy,r,S,L,J,T, Vsky0) * pow(h_width,2));	
			a[1] = 1. * (1. + 1./12. * v_HF(energy,r-h_width ,S,L,J,T, Vsky_) * pow(h_width,2));
			a[2] = 1. * (1. + 1./12. * v_HF(energy,r+h_width,S,L,J,T, Vsky1) * pow(h_width,2));
			return (a[0] * f0 - a[1] * f_) / a[2];
	}


	double v_HF(double energy, double r, double S, double L, double J, double T, double Vsky){
		if(T==-1./2.)
			return (energy-Vsky-centrifug_term(r,L)-potential_spin_orbit(r,S,L,J))/m_factor;
		else if(T==1./2.)
			return (energy-Vsky-centrifug_term(r,L)-potential_spin_orbit(r,S,L,J))/m_factor;
	}


	double skyrme(double rho, double rho_q, double rho_p, double rho_n){
		return  rho * ((t0/2.)*(2.+x0) + (2.+_a)*(t3/24.)*(2.+x3)*pow(rho,_a)) +
				rho_q * ((-t0/2.)*(2.*x0+1.) - (t3/12.)*(2.*x3 + 1.)*pow(rho,_a)) +
				_a * pow(rho,_a-1)*(-t3/24.)*(2.*x3 + 1.) * (pow(rho_p,2) + pow(rho_n,2));
	}

/*
	double normalise_HF(double eigen, double n_step_width_box, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1){
	vector<double> wfWork;
	wfWork.push_back(0.05);
	wfWork.push_back(0.1);

	double sum = 0.;
	for(int i=1; i<n_step_width_box+1; i++) wfWork.push_back(numerov_algorithm_HF(eigen, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, T, Vsky_, Vsky0, Vsky1));
	for(int i=1; i<n_step_width_box+1; i++) sum += h_width*pow(wfWork[i],2);

	wfWork.clear();
	return sum;
	}*/
	

/////////////////////OLD THINGS///////////////////////////


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
			return Z*e/(2.*R)*(3-pow(r/R,2));
		else
			return Z*e/r;
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

/*
	double normalise(double eigen, double n_step_width_box, double S, double L, double J, float T){
	vector<double> wfWork;
	wfWork.push_back(0.05);
	wfWork.push_back(0.1);

	double sum = 0.;
	for(int i=1; i<n_step_width_box+1; i++) wfWork.push_back(numerov_algorithm_woods(eigen, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, T));
	for(int i=1; i<n_step_width_box+1; i++) sum += h_width*pow(wfWork[i],2);

	wfWork.clear();
	return sum;
	}
	*/



