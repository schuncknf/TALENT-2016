#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <armadillo>
	#include <vector>

	#define Edown 0.
	#define Eup 10.
	#define width 35.
	#define h_width 0.01
	#define h_energy 0.0001

	#define prec 1E-10
	#define m_factor 20.75

	using namespace std;
	using namespace arma;


	double numerov_algorithm(double energy, double f_, double f0);

#endif
