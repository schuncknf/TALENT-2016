#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <vector>
	#include <fstream>

	#define selectFunc 4
	#define nProton 82
	#define nNeutron 126

	#define width_box 34.
	#define width 4.
	#define a_comp 0.67
	#define h_width 0.01
	#define V0 40.

	#define prec 1E-10
	#define m_factor 20.75

	using namespace std;

	double numerov_algorithm(double energy, double f_, double f0, double x);
	double V(double x);
	double infSW();
	double potential_well(double x);
	double woodsSaxon(double x);
	double V_so(double r, double l, double s);

#endif
