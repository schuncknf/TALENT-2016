#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <armadillo>
	#include <vector>

/*
	#include "TH1F.h"
	#include "TH2F.h"
	#include "TApplication.h"
	#include "TCanvas.h"
	#include "TGraph.h"
	#include "TAxis.h"
*/

	#define Edown 0.
	#define Eup 100.
	#define width 6.
	#define h_width 0.01
	#define h_energy 0.0005
	#define V0 40.

	#define width_box 34.

	#define prec 1E-10
	#define m_factor 20.75

	using namespace std;
	using namespace arma;


	double numerov_algorithm(double energy, double f_, double f0);
	double numerov_algorithm_finitewell(double energy, double f0, double f_,double x);
	double numerov_algorithm_woods(double energy, double f0, double f_,double x);
	double potential_well(double x);
	double potential_woods(double x);

#endif
