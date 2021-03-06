#include <iostream>	// cout cin
#include <stdio.h>	// sprintf
#include <string>	// strings
#include <fstream>	// input/output files
#include <sstream>	// ...
#include <iomanip>	// setprecision()
#include <stdlib.h>	// abs
#include <math.h>	// sin
#include <vector>	// vectors
#include <algorithm>	// min_element
#include <functional>	// 
#include <numeric>	// accumulate
#include "spherical.h"	// header file

using namespace std;	// no need to put std:: before any cout, etc.


struct state{
    std::vector<double> EFN;
    double n;
    double l;
    double j;
    double eEnergy;
};

struct compare{
    bool operator() (const state struct1, const state struct2)
    {
        return (struct1.eEnergy < struct2.eEnergy);
    }
};

// Main code
int main()
{   
    
    // Vectors of 'state' structs for proton and neutron states
    vector<state> protonStates;
    vector<state> neutronStates;
    
	// Open input file
	string line;
	ifstream ipFile;
	ipFile.open("input.dat");

	// Read in variables from file and store in vecot unless line is commented (starting with '#')
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
	//    selecctFunc = variable.at(0);	// Do not need, inherited from 1D code!
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

	int length = wBox/h;

	vector<double> totPot(length,0.0);
	vector<double> totProton(length,0.0);
	vector<double> totNeutron(length,0.0);

	vector<double> density(length,0.0);

	vector<double> engProton;
	vector<double> engNeutron;

	int nProt=0, nNeut=0;

	// Nice stuff for terminal output
// 	cout << "\n*************************************************************************" << endl;

	int nSteps = (eMax-eMin)/eStep;

	ofstream opFile2;
	opFile2.open("energies.dat");

	for(int isoSpin=0; isoSpin<2; isoSpin++)
	{
// 		if (isoSpin==0) cout << "*\tProtons \t\t\t\t\t\t\t*" << endl;
// 		else cout << "*\t\t\t\t\t\t\t\t\t*\n*\tNeutrons \t\t\t\t\t\t\t*" << endl;
// 		cout << "*\tEigen E [MeV]\tOrbital\t\tFilename\t\t\t*" << endl;

		if(isoSpin==0) opFile2 << "Protons:" << endl;
		if(isoSpin==1) opFile2 << "\nNeutrons:" << endl;

		for(int L=0; L<8; L++) // Loop up to L=7 ( i.e. j_15/2) orbitals
		{
			int orbCheck;
			if(L==0) orbCheck =1;
			else orbCheck=2;
			for(int spin=0; spin<orbCheck; spin++)
			{
				//eMin = variable.at(5);
				eMin = variable[5];
				double Spin;
				if(spin==0) Spin = -0.5;
				else Spin = 0.5;

				int n = 1;
                
                vector<double> potentialPrescan;
                for (int kk=1;kk<length;++kk){
                    potentialPrescan.push_back(V(kk*h,isoSpin,L,spin));
                }
                eMin = *min_element(potentialPrescan.begin(),potentialPrescan.end());
                nSteps = (eMax-eMin)/eStep;
//                 cout << eMin << endl;
                
				// Loop for the calculation
				for(int jj=0; jj<nSteps; jj++) // (eMax+eMin)
				{

					wf_val.push_back(wfStep1);
					wf_val.push_back(wfStep2);

					// Loop for calculating wavefunction values across the mesh, using the Numerov algorithm
					for(int i=2; i<wBox/h; i++)
					{
						//wf_val.push_back(numerovAlgorithm(eMin, wf_val.at(i+1), wf_val.at(i), i*h, isoSpin, L, spin));
						wf_val.push_back(numerovAlgorithm(eMin, wf_val[i-1], wf_val[i-2], i*h, isoSpin, L, spin));
					}

					// Assign value to wfLast
					//wfLast = wf_val.at(-1+wf_val.size());
					wfLast = wf_val[-1+wf_val.size()];

					// This is where will converge the calculation
					// Need if jj>0 condition otherwise will have no previous value for wavefunction
					if(jj>0)
					{
						//wfTmp = wf_val.at(-1+wf_val.size());
						wfTmp = wf_val[-1+wf_val.size()];

						// Check if there is a change in sign of wavefunction values at the limit of the box for
						// the most recent two different energies (wfPrev/wfLast negative if this is the case)
						if(wfPrev/wfLast < 0)
						{
							double eigenEng;	// initiate variable for energy eigenvalue

							// Check whether wavefunction at limit has pos. or neg. gradient at limit of box,
							// then calls the convergence function to calculate an energy eigenvalue
							if(wfPrev<wfLast)
							{
									eigenEng = converge(eMin-eStep,eMin,wfPrev,wfLast,isoSpin,L,spin);
							}
							else
							{
									eigenEng = converge(eMin,eMin-eStep,wfLast,wfPrev,isoSpin,L,spin);
							}
		
							// Clear previous wafefunction vector and initiate first two steps
							wf_val.clear();
							wf_val.push_back(wfStep1);
							wf_val.push_back(wfStep2);

							// Calculate eigenfunction using Numerov
							//for(int i=0; i<wBox/h; i++) wf_val.push_back(numerovAlgorithm(eigenEng, wf_val.at(i+1), wf_val.at(i), i*h, isoSpin, L, spin));
							for(int i=2; i<wBox/h; i++) wf_val.push_back(numerovAlgorithm(eigenEng, wf_val[i-1], wf_val[i-2], i*h, isoSpin, L, spin));

							// Calculate normalisation factor
							normFac = normalise(eigenEng, isoSpin, L, spin);

							double j=Spin+L;

							char orbital[16];
							if(L==0) sprintf(orbital,"%is_%1.0f.2",n,fabs(2*j));
							if(L==1) sprintf(orbital,"%ip_%1.0f.2",n,2*j);
							if(L==2) sprintf(orbital,"%id_%1.0f.2",n,2*j);
							if(L==3) sprintf(orbital,"%if_%1.0f.2",n,2*j);
							if(L==4) sprintf(orbital,"%ig_%1.0f.2",n,2*j);
							if(L==5) sprintf(orbital,"%ih_%1.0f.2",n,2*j);
							if(L==6) sprintf(orbital,"%ii_%1.0f.2",n,2*j);
							if(L==7) sprintf(orbital,"%ij_%1.0f.2",n,2*j);

							// Write results to file with filename including energy eigenvalue
							char file[512], filename[512], eigEng[512];

							if(isoSpin==0)
							{
								sprintf(filename,"results/proton_%s_%1.4f.dat",orbital,eigenEng);
								sprintf(file,"proton_%s_%1.4f.dat",orbital,eigenEng);
							}
							else
							{
								sprintf(filename,"results/neutron_%s_%1.4f.dat",orbital,eigenEng);
								sprintf(file,"neutron_%s_%1.4f.dat",orbital,eigenEng);
							}

							sprintf(eigEng,"%1.4f",eigenEng);

							ofstream opFile;
							opFile.open(filename);
							for(int k=0; k<wBox/h; k++) opFile << h*k << "\t" << sqrt(normFac)*wf_val.at(k) << "\t" << V(k*h, isoSpin, L, spin) << "\t" << woodsSaxon(k*h) << "\t" << centrifugal(k*h,L) << "\t" << spinOrbit(k*h,L,0.5) << "\t" << coulomb(k*h) << endl;
// 							for(int k=1; k<wBox/h; k++) opFile << h*k << "\t" << sqrt(normFac)*wf_val[k] << "\t" << pow((sqrt(normFac)*wf_val[k]/(k*h)),2) << endl;
							opFile.close();
                            
                            // Store state in 'protonStates' or 'neutronStates'
                            if(isoSpin==0){
                                state newProtonState;
                                
                                for(int k=0; k<wBox/h; k++) wf_val[k] = sqrt(normFac)*wf_val[k];
                                newProtonState.EFN = wf_val;
                                newProtonState.n = n;
                                newProtonState.l = L;
                                newProtonState.j = fabs(j);
                                newProtonState.eEnergy = eigenEng;
                                
                                protonStates.push_back(newProtonState);
                            }
                            else{
                                state newNeutronState;
                                
                                for(int k=0; k<wBox/h; k++) wf_val[k] = sqrt(normFac)*wf_val[k];
                                newNeutronState.EFN = wf_val;
                                newNeutronState.n = n;
                                newNeutronState.l = L;                                
                                newNeutronState.j = fabs(j);
                                newNeutronState.eEnergy = eigenEng;
                                
                                neutronStates.push_back(newNeutronState);
                            }
                                

							// Output energy eigenvalue and the name of the file results are saved to
							//cout << "*\t " << eigEng << "\t" << orbital << "\t\t" << file << "\t*"  << endl;
							n++;

							//opFile2 << fixed << setprecision(15) << orbital << "\t\t" << eigenEng << "\t" << (2*j)+1 << endl;
                            
							//vector<double> density(length,0.0);


						}
					}

					// Store the last value of the Numerov calculation 
					// so as to check whether change in sign in next iteration
					wfPrev = wfTmp;

					wf_val.clear(); // clear wf_val vector for next calculation with different trial energy
					eMin += eStep;	// increase energy.....will get rid of this
		
				}// Close loop over Energy mesh
				wfTmp =0;
				wfPrev=0;
			}	 // Close loop over +/- 1/2 spin
		
		}		 // Close looping over L values
	}			 // Close looping over isospin


// 	cout << "*************************************************************************\n" << endl;	// nice stuff for terminal output

	
    for(int isoSpin=0; isoSpin<2; isoSpin++)
    {
        for(unsigned int i=0; i<protonStates.size(); i++)
        {
            if(isoSpin==0 && nProt<=nProton)
            {
                nProt+= (2*protonStates[i].j)+1;
                for(int ii=1; ii<wBox/h; ii++)
                {
                    density[ii]   += ((2*protonStates[i].j)+ 1)*pow(protonStates[i].EFN[ii],2)/(4*M_PI*pow(ii*h,2));
                    totProton[ii] += ((2*protonStates[i].j)+ 1)*pow(protonStates[i].EFN[ii],2)/(4*M_PI*pow(ii*h,2));
                }
            }
        }
        
        for(unsigned int i=0; i<neutronStates.size(); i++)
        {      
            if(isoSpin==1 && nNeut<=nNeutron)
            {
                nNeut+= (2*neutronStates[i].j)+1;
                for(int ii=1; ii<wBox/h; ii++)
                {
                    density[ii]    += ((2*neutronStates[i].j)+ 1)*pow(neutronStates[i].EFN[ii],2)/(4*M_PI*pow(ii*h,2));
                    totNeutron[ii] += ((2*neutronStates[i].j)+ 1)*pow(neutronStates[i].EFN[ii],2)/(4*M_PI*pow(ii*h,2));
                }
            }
        }
    }
    
	cout << "\n*************************************************************************" << endl;
	cout << "*\tParticle totals\t\t\t\t\t\t\t*" << endl;
	cout << "*\tProtons:\t" << nProt << "\t\t\t\t\t\t*" << endl;
	cout << "*\tNeutrons:\t" << nNeut << "\t\t\t\t\t\t*" << endl;
    cout << "*\t\t\t\t\t\t\t\t\t*" << endl;
	
	ofstream opFile;
	opFile.open("densities.dat");
//	for(int i=0; i<wBox/h; i++) opFile << h*i << "\t" << totPot.at(i) << "\t" << totProton.at(i) << "\t" << totNeutron.at(i) << endl;
	for(int i=1; i<wBox/h; i++) opFile << h*i << "\t" << density[i] << "\t" << totProton[i] << "\t" << totNeutron[i] << endl;
	opFile.close();
    
    cout << "*************************************************************************" << endl;
    cout << "*\tParticle totals from density calculations" << "\t\t\t*" << endl;
	cout << "*\tTotal nucleons:\t" << totalMatterDensity(density) << "\t\t\t\t\t\t*" << endl;
	cout << "*\tTotal protons :\t" << totalMatterDensity(totProton) << " \t\t\t\t\t\t*" << endl;
	cout << "*\tTotal neutrons:\t" << totalMatterDensity(totNeutron) << " \t\t\t\t\t\t*" << endl;
    cout << "*\t\t\t\t\t\t\t\t\t*" << endl;
    
	density.clear();
	totProton.clear();
	totNeutron.clear();
	opFile2.close();
    
    sort(protonStates.begin(), protonStates.end(), compare());
    sort(neutronStates.begin(), neutronStates.end(), compare());
    
    string specNames = "spdfghij";
    
    cout << "*************************************************************************" << endl;
    cout << "*\tProton orbitals \t\t\t\t\t\t*" << endl;
    cout << "*\tEigen E [MeV]\tOrbital\t\tFilename\t\t\t*" << endl;
    cout << "*************************************************************************" << endl;
    for(unsigned int i=0;i<protonStates.size();++i){
        cout << "*\t" << setw(12) << setprecision(8) << fixed << protonStates[i].eEnergy << "\t" << setprecision(0) << neutronStates[i].n << specNames[protonStates[i].l] << " " << protonStates[i].j*2 << "/2\t\t\t\t\t\t*" << endl;
    }
    cout << "*\t\t\t\t\t\t\t\t\t*" << endl;
    
    cout << "*************************************************************************" << endl;
    cout << "*\tNeutron orbitals\t\t\t\t\t\t*" << endl;
    cout << "*\tEigen E [MeV]\tOrbital\t\tFilename\t\t\t*" << endl;
    cout << "*************************************************************************" << endl;
    for(unsigned int i=0;i<neutronStates.size();++i){
        cout << "*\t" << setw(12) << setprecision(8) << fixed << neutronStates[i].eEnergy << "\t" << setprecision(0) << neutronStates[i].n << specNames[neutronStates[i].l] << " " << neutronStates[i].j*2 << "/2\t\t\t\t\t\t*" << endl;
    }
    cout << "*\t\t\t\t\t\t\t\t\t*" << endl;
    cout << "*************************************************************************\n" << endl;
}


// Numerov algorithm function.
//
// INPUT ARGUMENTS:
//     E = trial energy
//   f_x = Wavefunction value from previous step ( f(x)   from TALENT school notes)
// f_x_h = Wavefunction value from two steps ago ( f(x-h) from TALENT school notes)
double numerovAlgorithm(double E, double f_x, double f_x_h, double r, int isoSpin, int L, int spin)
{
	double a[3];		// Array for Numerov coefficients
	double Vrplus = (E-V(r+h, isoSpin, L, spin)) / hm_fac;	// Numerov potential
	double Vrminus = (E-V(r-h, isoSpin, L, spin)) / hm_fac;	// Numerov potential
	double Vr0 = (E-V(r, isoSpin, L, spin)) / hm_fac;	// Numerov potential

//     if(r==0){
//         a[0] = 0.;	// Coeff. for f(x)
//     }
//     else{
//         a[0] = 2. * (1. - 5./12. * Vr0 * h * h);	// Coeff. for f(x)
//     }
//     if(r==0){
//         a[1] = 0.;	// Coeff. for f(x-h)
//     }
//     else if(r==h){
//         a[1] = 0.;	// Coeff. for f(x-h)
//     }
//     else{
//         a[1] = 1. * (1. + 1./12. * Vrminus * h * h);	// Coeff. for f(x-h)
//     }
//     a[2] = 1. * (1. + 1./12. * Vrplus * h * h);	// Coeff. for f(x+h)

    a[0] = 2. * (1. - 5./12. * Vr0 * h * h);
    a[1] = 1. * (1. + 1./12. * Vrminus * h * h);
    a[2] = 1. * (1. + 1./12. * Vrplus * h * h);

//     a[0] = 2. * (1. - 5./12. * Vrplus * h * h);
//     a[1] = 1. * (1. + 1./12. * Vrplus * h * h);
//     a[2] = 1. * (1. + 1./12. * Vrplus * h * h);
    
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
double converge(double eLo, double eHi, double wfLo, double wfHi, int isoSpin, int L, int spin)
{
	// Initiate variables that will be used to decide convergence
	long double trialE = (eHi+eLo)/2;
	vector<long double> wf_val;
	long double wfBarr=0, wfEnd = 10;

	long double diffE=eStep;
	// While condition to carry on bisection convergence energy difference between
	// trial energy and previous energies is small enough
	while (fabs(wfEnd)>1e-5 && diffE>convEng)
	{
		trialE = (eHi+eLo)/2;		// Trial energy based on average of two previous energies
		wf_val.push_back(wfStep1);	// Initiate Step 0 wavefunction value
		wf_val.push_back(wfStep2);	// Initiate Step 1 wavefunction value

		// Loop for calculating wavefunction values across the mesh, using the Numerov algorithm
		for(int i=2; i<wBox/h; i++)
		{
			//wf_val.push_back(numerovAlgorithm(trialE, wf_val.at(i+1), wf_val.at(i), i*h, isoSpin, L, spin));
			wf_val.push_back(numerovAlgorithm(trialE, wf_val[i-1], wf_val[i-2], i*h, isoSpin, L, spin));
		}

		// Work out which wavefunction is......I think this is where the convergence function is messing up.
		//if ( wf_val.at(-1+wf_val.size()) > 0)
		if ( wf_val[-1+wf_val.size()] > 0)
		{
			 eHi = trialE;
			//wfHi = wf_val.at(-1+wf_val.size());
			wfHi = wf_val[-1+wf_val.size()];
		}
		else
		{
			 eLo = trialE;
			//wfLo = wf_val.at(-1+wf_val.size());
			wfLo = wf_val[-1+wf_val.size()];
		}
		//wfEnd = wf_val.at(-1+wf_val.size());
		wfEnd = wf_val[-1+wf_val.size()];

		//wfBarr = wf_val.at((wBox+wWell)/(2*h));
		wfBarr = wf_val[(wBox+wWell)/(2*h)];

		wf_val.clear();	// clear wavefunction vector
		diffE = fabs(eHi-eLo);
	}

	if ( wfBarr/wfLo < 0 ) return eLo;
	else return eHi;
}

// Find factor for normalising wavefunction
double normalise(double eigenEng, int isoSpin, int L, int spin)
{
	vector<double> wfWork;
	wfWork.push_back(wfStep1);
	wfWork.push_back(wfStep2);

	double sum = 0;

	//for(int i=0; i<wBox/h; i++) wfWork.push_back(numerovAlgorithm(eigenEng, wfWork.at(i+1), wfWork.at(i), i*h, isoSpin, L, spin));
	//for(int i=0; i<(wBox/h)-1; i++) sum += h*pow(wfWork.at(i),2);

	for(int i=2; i<wBox/h; i++) wfWork.push_back(numerovAlgorithm(eigenEng, wfWork[i-1], wfWork[i-2], i*h, isoSpin, L, spin));

	for(int i=0; i<(wBox/h); i++) sum += h*pow((wfWork[i]),2);

	double normFac = 1/(sum);
	wfWork.clear();

	return normFac;
}

// Total potential
double V(double r, int isoSpin, int L, int spin)
{
//     if(r<r+h) r=h;
	if (isoSpin == 0)
		return woodsSaxon(r) + centrifugal(r,L) + spinOrbit(r,L,spin) + coulomb(r);
	else
		return woodsSaxon(r) + centrifugal(r,L) + spinOrbit(r,L,spin);
}

// Woods-Saxon potential
double woodsSaxon(double r)
{
	double WS;
	double         A = nProton + nNeutron;
	double      R = r0 * pow(A,1./3);
	double coeffV = -51. + 33.*((nNeutron-nProton)/A);

	WS = coeffV * (1 / ( 1 + exp( (r-R)/a0 )));

	return WS;
}

// Spin-orbit contribution
double spinOrbit(double r, int l, double spin)
{
	double s;
	if(l==0) return 0;
	else
	{
		if(spin==0) s = -0.5;
		else        s = 0.5;
		double    Vso = 0;
		int         A = nProton + nNeutron;
		double      R = r0 * pow(A,1./3);
		double      j = l+s;
		double coeffV = -51. + 33.*((nNeutron-nProton)/A);
		double   V_ls = Vls * 0.5 * coeffV * ( j*(j+1) - l*(l+1) - 3./4 ); // With underscore is total, without is strength coeff.

		Vso = V_ls * r0 * r0 * (1/r) * ( (exp((r+R)/a0)) / (a0*pow((exp(R/a0)+exp(r/a0)),2)) );

//		Vso = V_ls * r0 * r0 * (1/r) * ( (exp((r-R)/a0)) / (a0*pow((1+exp((r-R)/a0)),2)) );

		return Vso;
	}
}

// Centrifugal term
double centrifugal(double r, double l)
{
	double Vc = (hm_fac*l*(l+1)) / (r*r);
	return Vc;
}

// Coulomb potential
double coulomb(double r)
{
	double V_c;
	double Rp = r0 * pow(nProton,1./3);
	if(r<=Rp) V_c = (nProton * q / (2*Rp)) * (3 - pow(r/Rp,2));
	else V_c = nProton * q / r;

	return V_c;
}

double totalMatterDensity(vector<double> density)
{
	// Integrate total radial density 'density' to calculate A
	double A = 0.0;

	for(int i=1; i<wBox/h; i++)
	{
//		A += density.at(i);
		A += density[i]*4.*M_PI*(pow(i*h,2));
	}

	return h*A;
}