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
	#define _a 1
	#define t0 -1132.400
	#define t3 23610.40
	#define x0 0
	#define x3 0

	#define width_box 22.
	#define h_width 0.01
	#define r0 1.27
	#define e 1.439978

	#define A 208.
	#define N 126.
	#define Z 82.

	#define prec 1E-12
	#define m_factor 20.73553

	using namespace std;

	double numerov_algorithm_HF(double energy, double f0, double f_,double r, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1);
	double v_HF(double energy, double r, double S, double L, double J);
	double normalise_HF(double eigen, double n_step_width_box, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1);
	double centrifug_term(double r, double L);
	double potential_woods(double r);
	double potential_spin_orbit(double r, double S, double L, double J);
	double potential_coulomb(double r);
	double normalise(double eigen, double n_step_width_box, double S, double L, double J, float T);
	double skyrme(double rho, double rho_q, double rho_p, double rho_n);


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
