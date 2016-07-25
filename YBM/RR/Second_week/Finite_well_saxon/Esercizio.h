#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <armadillo>
	#include <vector>
	#include <fstream>
/*
	#include "TH1F.h"
	#include "TH2F.h"
	#include "TApplication.h"
	#include "TCanvas.h"
	#include "TGraph.h"
	#include "TAxis.h"
*/
	#define width_box 34.
	#define width 6.
	#define h_width 0.01
	#define V0 40.
	#define r0 1.27

	#define A 208
	#define N 82
	#define Z 127

	#define prec 1E-10
	#define m_factor 20.75

	using namespace std;

	double numerov_algorithm_woods_spin(double energy, double f0, double f_,double r,double S, double L);
	double potential_woods(double r);
	double potential_spin_orbit(double r, double S, double L);

#endif
