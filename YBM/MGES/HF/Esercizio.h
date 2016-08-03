#ifndef _Esercizio_H_	// Recursion guard

	#include <iostream>
	#include <cmath>
	#include <vector>
	#include <fstream>
	#include <ostream>
	#define _Esercizio_H_
	#include <sstream>
	#include <iomanip>
	#include <algorithm>

	using namespace std;

	// 'GLOBAL' VARIABLES //

	// Skyrme parameters for t0, t3 system
	#define _a 1
	#define t0 -1132.400
	#define t3 23610.40
	#define x0 0
	#define x3 0

	// System parameters
	#define width_box 22.
	#define h_width 0.01 // mesh spacing
	#define r0 1.27

	// Nucleons numbers
	#define A 16.
	#define N 8.
	#define Z 8.

	// Convergence parameters
	#define prec 1E-12	// precision
	#define minScanE -50. 	// energy above which to look for eigenfunctions
	#define maxProtScanE 10.	// energy below which to look for e'functions
	#define maxNeutScanE 0.

	#define integralPrec 1e-1 // iteration energy convergence

	// Physical constants
	#define m_factor 20.73553 // hBar^2/2m
	#define e 1.439978 // e^2

	// Function declarations
	double numerov_algorithm_HF(double energy, double f0, double f_,double r, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1);
	double v_HF(double energy, double r, double S, double L, double J, double T, double Vsky);
	double normalise_HF(double eigen, double n_step_width_box, double S, double L, double J, float T, double Vsky_, double Vsky0, double Vsky1);
	double centrifug_term(double r, double L);
	double potential_woods(double r);
	double potential_spin_orbit(double r, double S, double L, double J);
	double potential_coulomb(double r);
	double normalise(double eigen, double n_step_width_box, double S, double L, double J, float T);
	double skyrme(double rho, double rho_q, double rho_p, double rho_n);
	//int star();

	// Struct for storing all information about an e'function once found (and e'function itself)
	struct state{
		vector<double> wavefn;
    	double n, l, j;
    	double eig;
	};

	// Function to compare 2 eigenenergies from a 'state' object, for use in 'sort()'
	struct compare{
   		bool operator() (const state struct1, const state struct2){
   	    	return (struct1.eig < struct2.eig);
   		}
   	};


#endif
