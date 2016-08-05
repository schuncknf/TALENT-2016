#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <vector>
	#include <fstream>
	#include <ostream>
	#include <iomanip>
	#include <sstream>
	#include <algorithm>

	// Skyrme parameters for t0, t3 system
	#define _a 1.
	#define t0 -1132.400
	#define t3 23610.40
	#define x0 0.
	#define x3 0.
	#define x1 -0.317
	#define x2 -1.
	#define t1 484.23
	#define t2 -556.69

	//Costants and box features
	#define width_box 20.
	#define h_width 0.02
	#define r0 1.27
	#define e 1.439978
	#define m_factor 20.73553
	#define first_point 0.000005
	#define second_point 0.00001

	//Nucleon proprieties
	#define A 16.
	#define N 8.
	#define Z 8.
	#define max_n 2
	#define max_l 2

	//Convergence constants
	#define prec 1E-12
	#define HF_prec 1E-8

	using namespace std;


//Woods Saxon potential
	double numerov_algorithm_woods(double energy, double f0, double f_,double r, double S, double L, double J, float T);
	double v_neutron(double energy, double r,double S, double L, double J);
	double v_proton(double energy, double r,double S, double L, double J);
	double normalise(double eigen, double n_step_width_box, double S, double L, double J, float T);
	double centrifug_term(double r, double L);
	double potential_woods(double r);
	double potential_spin_orbit(double r, double S, double L, double J);
	double potential_coulomb(double r);

//Hartree Fock functions
	double numerov_algorithm_HF(double energy, double f0, double f_,double r, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1);
	double v_HF(double energy, double r, double S, double L, double J, double T, double Vsky);
	double skyrme(double rho, double rho_q, double rho_p, double rho_n, double rho1, double rho_q1, double rho2, double rho_q2, double r, double tau, double tau_q);


	struct state{
		vector<double> wavefn;
    	double n, l, j;
    	double eig;
	};
	
	struct compare{
   		bool operator() (const state struct1, const state struct2){
   	    	return (struct1.eig < struct2.eig);
   		}
   	};


#endif
