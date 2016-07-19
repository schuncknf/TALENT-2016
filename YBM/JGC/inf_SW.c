#include <iostream>	// cout cin
#include <math.h>	// sin
#include <vector>	// vectors
#include "inf_SW.h"	// header file

using namespace std;	// no need to put std:: before any cout, etc.

// Variables (NOTE: need to check keV and MeV being used correctly here....)
double     a = 35.;	// Box width  [fm]
double     l = 6.;	// Well width [fm] (don't use yet)
double     h = 0.01;	// Mesh width [fm]
double wellD = 1000000;	// Well depth [keV] (don't use yet)
double  eMin = 0.;	// Min limit  [keV]
double  eMax = 10.;	// Max limit  [keV]
double eStep = 0.0005;	// Step between different energies in brute force approach, needs getting rid of

double hm_fac = 20.75;	// units? MeV?

// Variable used to check if the wavefunction has crossed x axis
// between different energies
double wfPrev, wfLast;

// Main code
int main()
{
	double trialE = (eMax-eMin)/2;	// IS this the best way to start the bisection method?
	vector<double> wf_val;		// Vector for wave function values

	cout << "\n*************************************************************************" << endl;
	cout << "*\tE(i-1) [MeV]\tE(i) [MeV]\tWF(i-1) [a.u]\tWF(i) [a.u]\t*" << endl;

	// Loop for the calculation
	for(int jj=0; jj<(eMax-eMin)/eStep; jj++){
		
		wf_val.push_back(0);	// Initiate step 0 and step 1 WF values, if don't do this will not get past first calc.
		wf_val.push_back(0.5);	// NOTE: amplitude of calculated wavefunctions are arbitrary but defined by magnitude of
					//	 of this step, this is important when defining if calc. has converged.

		// Loop for calculating wavefunction values across the mesh, using the Numerov algorithm
		for(int i=0; i<a/h; i++){
			wf_val.push_back(numerovAlgorithm(eMin, wf_val.at(i+1), wf_val.at(i)));
		}
		
		// This is where will converge the calculation
		if(jj>0){
			// Assign value to wfLast
			wfLast = wf_val.at(-1+wf_val.size());

			// If there is a change in sign between different energies, cout info to terminal
			if( wfPrev/wfLast < 0 ) cout << "*\t" << eMin-eStep << "\t\t" << eMin << "\t\t" << wfPrev << "\t" << wfLast << "\t*" << endl;

/*
			// Below is where we can write results to file if converged......
			// still need to decide on what is "good enough"

			if( abs( (wf_val.at(-1+wf_val.size())) < 0.01)){
				char filename[512];
				sprintf(filename,"output_%1.3f.dat",eMin);
				ofstream opFile;
				opFile.open(filename);
				for(int k=0; k<a/h; k++) opFile << k << "\t" << wf_val.at(k) << endl;
				opFile.close();
			}
*/
		}

		// Store the last value of the Numerov calculation 
		// so as to check whether change in sign in next iteration
		wfPrev = wf_val.at(-1+wf_val.size());

		wf_val.clear(); // clear wf_val vector for next calculation with different trial energy
		eMin += eStep;	// increase energy.....will get rid of this
	}
	cout << "*************************************************************************\n" << endl;
}

// Numerov algorithm function.
//
// INPUT ARGUMENTS:
//     E = trial energy
//   f_x = Wavefunction value from previous step ( f(x) from TALENT school notes)
// f_x_h = Wavefunction value from two steps ago ( f(x) from TALENT school notes) 
double numerovAlgorithm(double E, double f_x, double f_x_h)
{
	double a[3];		// Array for Numerov coefficients
	double vx = E / hm_fac;	// Numerov potential

	a[0] = 2. * (1. - 5./12. * vx * h * h);	// Coeff. for f(x)
	a[1] = 1. * (1. + 1./12. * vx * h * h);	// Coeff. for f(x-h)
	a[2] = 1. * (1. + 1./12. * vx * h * h);	// Coeff. for f(x+h)

	double nuWF = ((a[0] * f_x) - (a[1] * f_x_h)) / a[2]; // Calculate wavefunction value
	return nuWF;	// Return nuWF value when function called
}
