#ifndef _Esercizio_H_
#define _Esercizio_H_

	#include <iostream>
	#include <cmath>
	#include <armadillo>
	#include <vector>

	#define Edown 0
	#define Eup 10
	#define prec 1E-10
	#define m_factor 20.75
	//#define h_bar 197.326947

	using namespace std;
	using namespace arma;

/*
	class Vettore{

	  public:
		Vettore();
		Vettore(int a);
		Vettore(Vettore &v);	
		~Vettore();
	
		void Setcomponent(int, double);
		void Insert();
		void Read();
		double Getcomponent(int n);

	  protected:
		int _n;
		double * _v;
	};

	Vettore wave_function(Vettore & x, double energy, double h);
*/

	double numerov_algorithm(double energy, double h, double f_, double f0);

#endif
