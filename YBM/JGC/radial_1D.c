#include <iostream>	// cout cin
#include <stdio.h>	// sprintf
#include <string>	// strings
#include <fstream>	// input/output files
#include <sstream>	// ...
#include <iomanip>	// setprecision()
#include <stdlib.h>	// abs
#include <math.h>	// sin
#include <vector>	// vectors
#include <numeric>	// accumulate
#include "radial_1D.h"	// header file

using namespace std;	// no need to put std:: before any cout, etc.

// Main code
int main()
{
	// Open input file
	string line;
	ifstream ipFile;
	ipFile.open("input.dat");

	// Read in data from file unless commented out (lines starting with '#' are comments
	while( getline( ipFile, line) )
	{
		if( line[0] != '#' )
		{
			istringstream iss(line);
			while(iss >> var)
			{
				variable.push_back(var);
			}
		}
	}

	// Set variables equal to those read in from input file
	selectFunc = variable.at(0);	// Selected potential
	      wBox = variable.at(1);	// Box width  [fm]
	     wWell = variable.at(2);	// Well width [fm] (don't use yet)
	         h = variable.at(3);	// Mesh width [fm]
	        V0 = variable.at(4);	// Well depth [MeV] (don't use yet)
	      eMin = variable.at(5);	// Min limit  [MeV]
	      eMax = variable.at(6);	// Max limit  [MeV]
	     eStep = variable.at(7);	// Step between different energies in brute force approach, needs getting rid of

	   wfStep1 = 0;			// First value for wavefunction, always set to 0
	   wfStep2 = variable.at(8);	// Second step along wavefunction
	   convEng = variable.at(9);	// Energy condition for convergence

	   nProton = variable.at(10);	// no. of protons for calculating Woods-Saxon
	  nNeutron = variable.at(11);	// no. of neutrons for calculating Woods-Saxon

	vector<double> wf_val;		// Vector for wave function values

	// Nice stuff for terminal output
	cout << "\n*************************************************" << endl;
	cout << "*\tEigen E [MeV]\tOutput filename\t\t*" << endl;

	int nSteps = (eMax-eMin)/eStep;

	// Loop for the calculation
	for(int jj=0; jj<nSteps; jj++) // (eMax+eMin)
	{
		
		wf_val.push_back(wfStep1);	// Initiate step 0 and step 1 WF values, if don't do this will not get past first calc.
		wf_val.push_back(wfStep2);	// NOTE: amplitude of calculated wavefunctions are arbitrary but defined by magnitude of
						//	 of this step, this is important when defining if calc. has converged.

		// Loop for calculating wavefunction values across the mesh, using the Numerov algorithm
		for(int i=0; i<wBox/h; i++)
		{
			wf_val.push_back(numerovAlgorithm(eMin, wf_val.at(i+1), wf_val.at(i), i*h));
		}

		// Assign value to wfLast
		wfLast = wf_val.at(-1+wf_val.size());

		// This is where will converge the calculation
		// Need if jj>0 condition otherwise will have no previous value for wavefunction
		if(jj>0)
		{
			wfTmp = wf_val.at(-1+wf_val.size());

			// Check if there is a change in sign of wavefunction values at the limit of the box for
			// the most recent two different energies (wfPrev/wfLast negative if this is the case)
			if(wfPrev/wfLast < 0)
			{
				double eigenEng;	// initiate variable for energy eigenvalue

				// Check whether wavefunction at limit has pos. or neg. gradient at limit of box,
				// then calls the convergence function to calculate an energy eigenvalue
				if(wfPrev<wfLast)
				{
						eigenEng = converge(eMin-eStep,eMin,wfPrev,wfLast);
				}
				else
				{
						eigenEng = converge(eMin,eMin-eStep,wfLast,wfPrev);
				}
			
				// Clear previous wafefunction vector and initiate first two steps
				wf_val.clear();
				wf_val.push_back(wfStep1);
				wf_val.push_back(wfStep2);

				// Calculate eigenfunction using Numerov
				for(int i=0; i<wBox/h; i++) 
				{
					if (i < (wBox+wWell)/(2*h)) wf_val.push_back(numerovAlgorithm(eigenEng, wf_val.at(i+1), wf_val.at(i), i*h));
					else if ( (numerovAlgorithm(eigenEng, wf_val.at(i+1), wf_val.at(i), i*h))/wf_val.at((wBox+wWell)/(2*h))>0) wf_val.push_back(numerovAlgorithm(eigenEng, wf_val.at(i+1), wf_val.at(i), i*h));
					else wf_val.push_back(0);
				}

				double normFac = normalise(eigenEng);

				// Write results to file with filename including energy eigenvalue
				char filename[512], eigEng[512];
				sprintf(filename,"output_%1.4f.dat",eigenEng);
				sprintf(eigEng,"%1.4f",eigenEng);
				ofstream opFile;
				opFile.open(filename);
				for(int k=0; k<wBox/h; k++) opFile << h*k << "\t" << sqrt(normFac)*wf_val.at(k) << endl;
				opFile.close();

				// Output energy eigenvalue and the name of the file results are saved to
				//cout << "*\t" << setprecision(10) << eigEng << "\t\t" << filename << "\t*"  << endl;
			}
		}

		// Store the last value of the Numerov calculation 
		// so as to check whether change in sign in next iteration
		wfPrev = wfTmp;

		wf_val.clear(); // clear wf_val vector for next calculation with different trial energy
		eMin += eStep;	// increase energy.....will get rid of this
	}
	cout << "*************************************************\n" << endl;	// nice stuff for terminal output
}

// Numerov algorithm function.
//
// INPUT ARGUMENTS:
//     E = trial energy
//   f_x = Wavefunction value from previous step ( f(x)   from TALENT school notes)
// f_x_h = Wavefunction value from two steps ago ( f(x-h) from TALENT school notes)
double numerovAlgorithm(double E, double f_x, double f_x_h, double x)
{
	double a[3];		// Array for Numerov coefficients
	double vx = (E-V(x)) / hm_fac;	// Numerov potential

	a[0] = 2. * (1. - 5./12. * vx * h * h);	// Coeff. for f(x)
	a[1] = 1. * (1. + 1./12. * vx * h * h);	// Coeff. for f(x-h)
	a[2] = 1. * (1. + 1./12. * vx * h * h);	// Coeff. for f(x+h)

	double nuWF = ((a[0] * f_x) - (a[1] * f_x_h)) / a[2]; // Calculate wavefunction value

	return nuWF;	// Return nuWF value when function called
}

// Function to converge on energy eigenvalue
//
// INPUT ARGUMENTS:
//  eLo = the energy of the step with the lowest wavefunction value
//  eHi = the energy of the step with the highest wavefunction value
// wfLo = lower wavefunction value of the two energy steps
// wfHi = higher wavefunction value of the two energy steps
double converge(double eLo, double eHi, double wfLo, double wfHi)
{
	long double trialE = (eHi+eLo)/2;
	vector<long double> wf_val;
	long double wfBarr=0, wfEnd = 10;

//	cout << grad << "\t" << wfHi << "\t" << wfLo << endl;
	long double diffE=eStep;
	// While condition to carry on bisection convergence energy difference between
	// trial energy and previous energies is small enough
	while (fabs(wfEnd)>1e-7 && diffE>1e-14)
	{
		trialE = (eHi+eLo)/2;		// Trial energy based on average of two previous energies
		wf_val.push_back(wfStep1);	// Initiate Step 0 wavefunction value
		wf_val.push_back(wfStep2);	// Initiate Step 1 wavefunction value

		// Loop for calculating wavefunction values across the mesh, using the Numerov algorithm
		for(int i=0; i<wBox/h; i++)
		{
			wf_val.push_back(numerovAlgorithm(trialE, wf_val.at(i+1), wf_val.at(i), i*h));
		}

		// Work out which wavefunction is......I think this is where the convergence function is messing up.
		if ( wf_val.at(-1+wf_val.size()) > 0)
		{
			 eHi = trialE;
			wfHi = wf_val.at(-1+wf_val.size());
		}
		else
		{
			 eLo = trialE;
			wfLo = wf_val.at(-1+wf_val.size());
		}
		wfEnd = wf_val.at(-1+wf_val.size());

		wfBarr = wf_val.at((wBox+wWell)/(2*h));

		wf_val.clear();	// clear wavefunction vector
		diffE = fabs(eHi-eLo);
	}

	if ( wfBarr/wfLo < 0 ) return eLo;
	else return eHi;

//	return eHi;	// return the trial energy value that meets the tolerance condition
}

// Find factor for normalising wavefunction
double normalise(double eigenEng)
{
	vector<double> wfWork;
	wfWork.push_back(wfStep1);
	wfWork.push_back(wfStep2);

	double sum = 0;

//	for(int i=0; i<wBox/h; i++) wfWork.push_back(fabs((numerovAlgorithm(eigenEng, wfWork.at(i+1), wfWork.at(i), i*h))));

	for(int i=0; i<wBox/h; i++) 
	{
		if (i < (wBox+wWell)/(2*h)) wfWork.push_back(numerovAlgorithm(eigenEng, wfWork.at(i+1), wfWork.at(i), i*h));
		else if ( (numerovAlgorithm(eigenEng, wfWork.at(i+1), wfWork.at(i), i*h))/wfWork.at((wBox+wWell)/(2*h))>0) wfWork.push_back(numerovAlgorithm(eigenEng, wfWork.at(i+1), wfWork.at(i), i*h));
		else wfWork.push_back(0);
	}

	for(int i=0; i<(wBox/h)-1; i++) sum += h*pow(wfWork.at(i),2);

	double normFac = 1/sum;
	wfWork.clear();

	return normFac;
}

double V(double x)
{
	if (selectFunc == 1) return infSW();
	if (selectFunc == 2) return finSW(x);
	if (selectFunc == 3) return woodsSaxon(x);
}

double testPot(int t)
{
	ofstream tFile;
	tFile.open("test_pot.dat");

	if (t == 2) for(int i=0; i<wBox/h; i++) tFile << h*i << "\t" << finSW(h*i) << endl;
	if (t == 3) for(int i=0; i<wBox/h; i++) tFile << h*i << "\t" << woodsSaxon(h*i) << endl;

	tFile.close();
}

// Infinite square well potential
double infSW()
{
	return 0;
}

// Finite square well
double finSW(double x)
{
	if( (x<((wBox-wWell)/2)) || (x>((wBox+wWell)/2))  ) 
		return 0;
	else
		return -V0;
}

// Woods-Saxon potential
double woodsSaxon(double x)
{
	double WS;
	int         A = nProton + nNeutron;
	double     r0 = 1.27;
	double      R = r0 * pow(A,1./3);
	double a_comp = 0.67;
	double coeffV = -51 + 33*((nNeutron-nProton)/A);

	if(x>wBox/2)
	{
		WS = coeffV * (1 / ( 1 + exp( ((x-(wBox/2))-R)/a_comp )));
	}
	else 
	{
		WS = coeffV * (1 / ( 1 + exp( ((-x+(wBox/2))-R)/a_comp )));
	}
	return WS;
}
