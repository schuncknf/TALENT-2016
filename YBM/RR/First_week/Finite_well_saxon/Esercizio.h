#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <armadillo>
	#include <vector>
	#include <fstream>

//Box and well features
	#define width_box 34.
	#define width 6.
	#define h_width 0.01
	#define V0 40.

//Constants and precisions
	#define prec 1E-10
	#define m_factor 20.75

//Eigenvalues to plot, up to 2
	#define n 3


	using namespace std;

	double numerov_algorithm(double energy, double f_, double f0);
	double numerov_algorithm_finitewell(double energy, double f0, double f_,double x);
	double numerov_algorithm_woods(double energy, double f0, double f_,double x);
	double potential_well(double x);
	double potential_woods(double x);

#endif
