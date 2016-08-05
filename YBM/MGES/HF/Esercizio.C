#include "Esercizio.h"
#include "initial.h"
#include "star.h"

int main(){

	// Declaration of constants and variables
	int n_step_width_box = width_box/h_width;
	vector<double> wave_val;	// Vector for storing e'functions during
	vector<double> wfWork;	// Vector for storing e'function and then normalising
	double S = 1./2.;	// Spin
	vector<double> U_skyrme_p, U_skyrme_n, U_cou_p;	// Neutron and proton potential vectors

	double wavefn0 = 0.000005;
	double wavefn1 = 0.00001;

	// Variables for total system of energy at current and previous iteration
	// Initialise with values such that 'while' loop actually begins
   	double integral=0., old_integral=1e10;

	// Arrays for densities
	double * density_neutron = new double [n_step_width_box+1];
	double * density_proton = new double [n_step_width_box+1];
	double * K_skyrme = new double [n_step_width_box+1];
	double * E_skyrme = new double [n_step_width_box+1];

	// 'state' structs to store neutron and proton eigenstates
	state appoggio;
	vector<state> neutronstates;
	vector<state> protonstates;


	// Read initial densities from files

    ifstream myneutron, myproton;
    myneutron.open("Oxy_density_neutron.txt");
    myproton.open("Oxy_density_proton.txt");

	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
		density_neutron[i] = 0.;
	}

	double radius, fileDensity;
	// Neutrons
	int fileLine = 0;
    while(myneutron >> radius >> fileDensity){
		density_neutron[fileLine] = fileDensity;
		fileLine++;
	}
	myneutron.close();

	// Protons
    fileLine = 0;
    while(myproton >> radius >> fileDensity){
		density_proton[fileLine] = fileDensity;
		fileLine++;
	}
	myproton.close();
	// for(int i=0; i<n_step_width_box; i++) cout << density_proton[i] << "\t" << density_neutron[i] << endl;


	// Hartree Fock iteration loop
	int iteration = 0;
	while(fabs(old_integral - integral) > integralPrec){

		++iteration;
		cout << "\nIteration number =\t" << iteration << endl;

		neutronstates.clear();
		protonstates.clear();


		// Calculate Skyrme potentials from total, proton, and neutron densities
	    U_skyrme_p.clear();
		U_skyrme_n.clear();
		U_cou_p.clear();

		for(int i=0; i<n_step_width_box; i++){
			double r  = density_proton[i] + density_neutron[i];
			double rp = density_proton[i];
			double rn = density_neutron[i];
			// Inputs are rho, rho_q, rho_p, rho_n
			U_skyrme_p.push_back(skyrme(r,rp,rp,rn));
			U_skyrme_n.push_back(skyrme(r,rn,rp,rn));


			double integ1 = 0., integ2 = 0., integ3 = 0., exchange = 0.;

			for(int jj=0;jj<=i;++jj){
				integ1 += density_proton[jj]*pow(jj*h_width+h_width,2);
				integ2 += density_proton[jj]*(jj*h_width+h_width);
			}
			for(int jj=0;jj<=width_box/h_width;++jj){
				integ3 += density_proton[jj]*(jj*h_width+h_width);
			}
			exchange = pow(density_proton[i],1./3.);
			U_cou_p.push_back( h_width* (4.*M_PI*e* (integ1/(i*h_width+h_width) - integ2 + integ3) - exchange*e*pow(3./M_PI,1./3.)) );
		}

		// Neutron eigenstates
		double isospin = -1./2.;
		// Loop over number of nodes (i.e. quantum no "n")
		for(int node=0; node<3; node++){
			// Loop over angular momentum
			for(int L=0; L<8;L++){
				// Loop over spin up/down
				for(int h=1; h<3; h++){

					// Calculate j=l+s
					if(L==0) h=2;
					double J = abs(S+pow(-1,h)*L);

					int nodecount=0;
					double Etrial=0.;
					double Eup = 0.;
					double Edown = minScanE;

					// Convergence loop
					while(abs(Eup-Edown)>prec){

						// Integrate with trial energy
						Etrial = (Eup+Edown)/2.;
						wave_val.push_back(wavefn0);
						wave_val.push_back(wavefn1);
						for(int i=1; i<n_step_width_box; i++){
							wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, isospin, U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
						}

						// Count nodes in resulting wavefunction
						nodecount = 0;
						for(int i=1; i<n_step_width_box+2; i++){
							if(wave_val[i-1]*wave_val[i]<0.)
							nodecount++;
						}

						// Adjust energy window according to no. of nodes
						if(nodecount>node){
							Eup = Etrial;
						}
						else if(nodecount<node||nodecount==node){
							Edown = Etrial;
						}
						wave_val.clear();
					}

					// Store converged eigenstate, if e'energy in suitable range
					if(Etrial<maxNeutScanE && Etrial>minScanE){

						// Calculate normalisation factor 'norm'
						wfWork.push_back(wavefn0);
						wfWork.push_back(wavefn1);
						double norm = 0.;
						for(int i=1; i<n_step_width_box+1; i++){
							wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, isospin, U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
						}
						for(int i=1; i<n_step_width_box+1; i++) norm += h_width*pow(wfWork[i],2);
						wfWork.clear();

						// Store e'state
						appoggio.wavefn.push_back(wavefn0/sqrt(norm));
						appoggio.wavefn.push_back(wavefn1/sqrt(norm));
						for(int i=1; i<n_step_width_box; i++){
							appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, isospin, U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1])/sqrt(norm));
						}
						appoggio.eig = Etrial;
						appoggio.l = L;
						appoggio.j = J;
						neutronstates.push_back(appoggio);
						appoggio.wavefn.clear();
					}

				}	// Spin loop end
			}	// Angular momentum loop end
		}	// Nodes loop end



		// Proton eigenstates
		isospin = 1./2.;

		// Loop over number of nodes (i.e. quantum no "n")
		for(int node=0; node<3; node++){
			// Loop over angular momentum
			for(int L=0; L<8;L++){
				// Loop over spin up/down
				for(int h=1; h<3; h++){

					// Calculate j=l+s
					if(L==0) h=2;
					double J = abs(S+pow(-1,h)*L);

					int nodecount=0;
					double Etrial=0.;
					double Eup = 0.;
					double Edown = minScanE;

					// Convergence loop
					while(abs(Eup-Edown)>prec){

						// Integrate with trial energy
						Etrial = (Eup+Edown)/2.;
						wave_val.push_back(wavefn0);
						wave_val.push_back(wavefn1);
						for(int i=1; i<n_step_width_box; i++){
							wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_p[i-1]+U_cou_p[i-1], U_skyrme_p[i]+U_cou_p[i], U_skyrme_p[i+1]+U_cou_p[i+1]));
						}

						// Count nodes in resulting wavefunction
						nodecount = 0;
						for(int i=1; i<n_step_width_box+2; i++){
							if(wave_val[i-1]*wave_val[i]<0.)
							nodecount++;
						}

						// Adjust energy window according to no. of nodes
						if(nodecount>node){
							Eup = Etrial;
						}
						else if(nodecount<node||nodecount==node){
							Edown = Etrial;
						}
						wave_val.clear();
					}

					// Store converged eigenstate, if e'energy in suitable range
					if(Etrial<maxProtScanE && Etrial>minScanE){

						// Calculate normalisation factor 'norm'
						wfWork.push_back(wavefn0);
						wfWork.push_back(wavefn1);
						double norm = 0.;
						for(int i=1; i<n_step_width_box; i++) wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_p[i-1]+U_cou_p[i-1], U_skyrme_p[i]+U_cou_p[i], U_skyrme_p[i+1]+U_cou_p[i+1]));
						for(int i=1; i<n_step_width_box+1; i++) norm += h_width*pow(wfWork[i],2);
						wfWork.clear();

						// Store e'state
						appoggio.wavefn.push_back(wavefn0);
						appoggio.wavefn.push_back(wavefn1);
						for(int i=1; i<n_step_width_box; i++){
							appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_p[i-1]+U_cou_p[i-1], U_skyrme_p[i]+U_cou_p[i], U_skyrme_p[i+1]+U_cou_p[i+1])/sqrt(norm));
						}
						appoggio.eig = Etrial;
						appoggio.l = L;
						appoggio.j = J;
						protonstates.push_back(appoggio);
						appoggio.wavefn.clear();
					}

				}	// Spin loop end
			}	// Angular momentum loop end
		}	// Nodes loop end


		// Sort vectors of neutron and proton states by e'energies

		int n_neutron = 0, n_levels = 0;
		sort(neutronstates.begin(), neutronstates.end(), compare());
		cout << "J value\t\t\t" << "Eigenvalue" << endl;
		for(int i=0; i<neutronstates.size(); i++){
			if(n_neutron<N){
				cout  <<  neutronstates[i].j*2 << "/2\t\t\t" << setprecision(12) << neutronstates[i].eig  <<  endl;
				n_neutron += (2*neutronstates[i].j+1);
				n_levels ++;
			}
		}
		cout << "Il numero di neutroni richiesto è:" << n_neutron  <<  endl;

		int n_proton = 0, p_levels = 0;
		sort(protonstates.begin(), protonstates.end(), compare());
		cout << "J value\t\t\t" << "Eigenvalue" << endl;
		for(int i=0; i<protonstates.size(); i++){
			if(n_proton<Z){
				cout  <<  protonstates[i].j*2 << "/2\t\t\t" << setprecision(12) << protonstates[i].eig  <<  endl;
				n_proton += (2*protonstates[i].j+1);
				p_levels ++;
			}
		}
		cout<<"Il numero di protoni richiesto è:"<<n_proton<<endl;


		// Calculate new densities for neutrons and protons
		for(int i=0; i<n_step_width_box+1; i++){
			density_neutron[i] = 0.;
		}
		for(int k=0; k<n_levels; k++){
			for(int i=2; i<n_step_width_box; i++)
				density_neutron[i] += (2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(neutronstates[k].wavefn[i],2);
		}

		for(int i=0; i<n_step_width_box+1; i++){
			density_proton[i] = 0.;
		}
		for(int k=0; k<p_levels; k++){
			for(int i=2; i<n_step_width_box; i++)
				density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
		}


		// Calculate total energy of system using kinetic and Skyrme densities

		// Calculate the kinetic density
		for(int i=0; i<n_step_width_box+1; i++){
			K_skyrme[i] = 0.;
		}
		for(int k=0; k<n_levels; k++){
			for(int i=2; i<n_step_width_box; i++)
				K_skyrme[i] += m_factor*(2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*(pow((neutronstates[k].wavefn[i]-neutronstates[k].wavefn[i-1])/h_width-neutronstates[k].wavefn[i]/(i*h_width+h_width),2)+
										neutronstates[k].l*(neutronstates[k].l+1)/pow(i*h_width+h_width,2)*pow(neutronstates[k].wavefn[i],2) ) 	// kinetic energy of neutrons
							+ m_factor*(2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*(pow((protonstates[k].wavefn[i]-protonstates[k].wavefn[i-1])/h_width-protonstates[k].wavefn[i]/(i*h_width+h_width),2)+
										protonstates[k].l*(protonstates[k].l+1)/pow(i*h_width+h_width,2)*pow(protonstates[k].wavefn[i],2) ); 	// kinetic energy of protons
		}

		// Calculate the skyrme density
		for(int i=0; i<n_step_width_box+1; i++){
			E_skyrme[i] = 0.;
		}
		for(int i=2; i<n_step_width_box; i++){
			double r  = density_proton[i] + density_neutron[i];
			double rp = density_proton[i];
			double rn = density_neutron[i];
			E_skyrme[i] += 	1./2.*t0*((1+x0/2.)*pow(r,2)-(x0+1./2.)*(pow(rp,2)+pow(rn,2)))+
									1./12.*t3*pow(r,_a)*((1+x3/2.)*pow(r,2)-(x3+1./2.)*(pow(rp,2)+pow(rn,2)));
				//cout<<"K_skyrme:"<<"\t"<<K_skyrme[i]<<"\t"<<"E_skyrme:"<<"\t"<<E_skyrme[i]<<endl;
		}

		// Integral densities to get new total energy 'integral'
		old_integral = integral;
		double integralKin = 0., integralSky = 0.;
		for(int i=0; i<n_step_width_box+1; i++){
			integralKin += 4*M_PI*pow(i*h_width+h_width,2)*K_skyrme[i]*h_width;
			integralSky += 4*M_PI*pow(i*h_width+h_width,2)*E_skyrme[i]*h_width;
		}
		integral = integralKin + integralSky;

		cout << "Kinetic energy:\t" << integralKin << "\nSkyrme energy:\t" << integralSky << endl;
		cout << "Total energy:\t" << integral << "\nChange from previous iteration:\t" << old_integral-integral << endl;
	}



	int n_neutron = 0, n_levels = 0;
	int n_proton = 0, p_levels = 0;

	cout << "\n\nJ value\t\t\t" << "Eigenvalue" << endl;
	for(int i=0; i<neutronstates.size(); i++){
		if(n_neutron<N){
			cout  <<  neutronstates[i].j*2 << "/2\t\t\t" << setprecision(12) << neutronstates[i].eig  <<  endl;
			n_neutron += (2*neutronstates[i].j+1);
			n_levels ++;
		}
	}
	cout << "Il numero di neutroni richiesto è:" << n_neutron  <<  endl;

	cout << "J value\t\t\t" << "Eigenvalue" << endl;
	for(int i=0; i<protonstates.size(); i++){
		if(n_proton<Z){
			cout  <<  protonstates[i].j*2 << "/2\t\t\t" << setprecision(12) << protonstates[i].eig  <<  endl;
			n_proton += (2*protonstates[i].j+1);
			p_levels ++;
		}
	}
	cout<<"Il numero di protoni richiesto è:"<<n_proton<<endl;

	cout << "Total energy:\t" << integral << endl;

	//star();


	//Routine to obtain proton density, and the sum of both

	// double * density_proton = new double [n_step_width_box+1];
	// for(int i=0; i<n_step_width_box+1; i++){
	// 	density_proton[i] = 0.;
	// }
	// for(int k=0; k<p_levels-7; k++){
	// 	for(int i=2; i<n_step_width_box; i++)
	// 		density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	// }

	ofstream myfile_density_neutron, myfile_density_proton, myfile_density_np;
	myfile_density_neutron.open("density_neutron.txt");
	myfile_density_proton.open("density_proton.txt");
	myfile_density_np.open("density_np.txt");

	for(int i=2; i<n_step_width_box+1; i++){
		myfile_density_neutron << i*h_width +h_width<< "   " << density_neutron[i] << endl;
		myfile_density_proton << i*h_width +h_width<< "   " << density_proton[i] << endl;
		myfile_density_np << i*h_width +h_width << "   " << density_neutron[i] + density_proton[i] << endl;
	}
	myfile_density_neutron.close();
	myfile_density_proton.close();
	myfile_density_np.close();






//USEFUL ROUTINES

//This is just for plotting all wavefunction
	/*
cout<<n_levels<<endl;
	char orbital[26];
	ofstream myfile_wave;
	for(int k=0; k<n_levels; k++){
					sprintf(orbital,"wf_neutron%1.4f.dat",neutronstates[k].eig);
					myfile_wave.open(orbital);
		for(int i=2; i<n_step_width_box+1; i++){
					myfile_wave << i*h_width +h_width<< "   " << neutronstates[k].wavefn[i]/(i*h_width +h_width) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
				}
				myfile_wave.close();
			}*/



//This can be useful to plot a potential

	ofstream myfile;
	myfile.open("Eigen.txt");

	for(int i=1; i<n_step_width_box; i++){
		//wave_val.push_back(potential_woods(i*h_width +h_width)-centrifug_term(i*h_width +h_width, 0));
		myfile <<i*h_width +h_width<< "   " << U_skyrme_n[i]<< endl;
	//	myfile <<i*h_width +h_width<< "\t" << U_skyrme_p[i] << "\t" << U_skyrme_n[i] << endl;
	}

	myfile.close();



//Plot the proton and neutron central potential
/*
	ofstream myfileN, myfileP;
	myfileN.open("EigenN.txt");
	myfileP.open("EigenP.txt");

	vector<double> Proton_pot;
	vector<double> Neutron_pot;

	for(int i=1; i<n_step_width_box; i++){
		Proton_pot.push_back(potential_woods(i*h_width+h_width));
		Neutron_pot.push_back(potential_woods(i*h_width+h_width)+potential_coulomb(i*h_width+h_width));
		myfileN <<i*h_width +h_width<< "   " << Proton_pot[i] << endl;
		myfileP <<i*h_width +h_width<< "   " << Neutron_pot[i] << endl;
	}

	myfileN.close();
	myfileP.close();
*/

	return 0;
}
