#include"Esercizio.h"


	double numerov_algorithm_woods(double energy, double f0, double f_,double r, double S, double L, double J, float T){
		double v=0., a[3];
		if(T==-1./2.){
			a[0] = 2. * (1. - 5./12. * v_neutron(energy,r,S,L,J) * pow(h_width,2));	
			a[1] = 1. * (1. + 1./12. * v_neutron(energy,r-h_width ,S,L,J) * pow(h_width,2));
			a[2] = 1. * (1. + 1./12. * v_neutron(energy,r+h_width,S,L,J) * pow(h_width,2));
			return (a[0] * f0 - a[1] * f_) / a[2];
		}
		else if(T==1./2.){
			a[0] = 2. * (1. - 5./12. * v_proton(energy,r,S,L,J) * pow(h_width,2));
			a[1] = 1. * (1. + 1./12. * v_proton(energy,r-h_width,S,L,J) * pow(h_width,2));
			a[2] = 1. * (1. + 1./12. * v_proton(energy,r+h_width,S,L,J) * pow(h_width,2));
			return (a[0] * f0 - a[1] * f_) / a[2];
		}
	}

	double v_neutron(double energy, double r,double S, double L, double J){
		return (energy-potential_woods(r)-centrifug_term(r,L)-potential_spin_orbit(r,S,L,J))/m_factor;
	}

	double v_proton(double energy, double r,double S, double L, double J){
		return (energy-potential_woods(r)-centrifug_term(r,L)-potential_spin_orbit(r,S,L,J)-potential_coulomb(r))/m_factor;
	}

	double normalise(double eigen, double n_step_width_box, double S, double L, double J, float T){
	vector<double> wfWork;
	wfWork.push_back(first_point);
	wfWork.push_back(second_point);

	double sum = 0.;
	for(int i=1; i<n_step_width_box+1; i++) wfWork.push_back(numerov_algorithm_woods(eigen, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, T));
	for(int i=1; i<n_step_width_box+1; i++) sum += h_width*pow(wfWork[i],2);

	wfWork.clear();
	return sum;
}

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
			return (energy-Vsky-centrifug_term(r,L)/*-potential_spin_orbit(r,S,L,J)*/)/m_factor;
		else if(T==1./2.)
			return (energy-Vsky-centrifug_term(r,L)/*-potential_spin_orbit(r,S,L,J)*/)/m_factor;
	}


	double skyrme(double * density_q, double * density_n, double * density_p, int i, double r, double tau, double tau_q){
		double d_density = ((density_p[i+1]+density_n[i+1])-(density_p[i]+density_n[i]))/h_width;
		double d_density_q = (density_q[i+1]-density_q[i])/h_width;
		double d2_density = (2*(density_p[i+1]+density_n[i+1]) - 2*(density_p[i]+density_n[i]) - 2*d_density*h_width)/pow(h_width,2);
		double d2_density_q = (2*density_q[i+1] - 2*density_q[i] - 2*d_density_q*h_width)/pow(h_width,2);

		return  (density_n[i]+density_p[i]) * ((t0/2.)*(2.+x0) + (2.+_a)*(t3/24.)*(2.+x3)*pow(density_n[i]+density_p[i],_a)) +
				density_q[i] * ((-t0/2.)*(2.*x0+1.) - (t3/12.)*(2.*x3 + 1.)*pow(density_n[i]+density_p[i],_a)) +
				_a * pow(density_n[i]+density_p[i],_a-1)*(-t3/24.)*(2.*x3 + 1.) * (pow(density_p[i],2) + pow(density_n[i],2))+
				tau*1./8.*(t1*(2.+x1)+t2*(2.+x2))-
				1./16.*(3.*t1*(2.+x1)-t2*(2.+x2))*(d2_density+2.*d_density/r)-
				1./8.*(t1*(2.*x1+1.)-t2*(2.*x2+1.))*tau_q-
				1./8.*(t1*x1+t2*x2)*(d2_density_q+2.*d_density_q/r);
	}


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



