#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <vector>
	#include <fstream>

	#define width_box 34.
	#define width 6.
	#define h_width 0.001
	#define V0 40.

	#define prec 1E-10
	#define m_factor 20.75
    
    #define r0 1.27
    #define a_comp 0.67
    #define Vls 0.44
    
    #define nProton 82
    #define nNeutron 126

	using namespace std;

	double numerov_algorithm(double energy, double f_, double f0, double x, int potType, int l, double s);
    double potential(double x, int potType, int l, double s);

#endif
