#include "Esercizio.h"

int main(){


//Some costants and variables

	double Edown = 0.;
	double Eup = -50.;
	int n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	double S = 1./2.;
    vector<state> neutronstates;
    vector<double> U_skyrme_p, U_skyrme_n;
   	double integral=0., old_integral;

	double * density_neutron = new double [n_step_width_box+1];
	double * density_proton = new double [n_step_width_box+1];
	double * K_skyrme = new double [n_step_width_box+1];
	double * E_skyrme = new double [n_step_width_box+1];

//I get the first density from files

    ifstream myneutron, myproton;
    myneutron.open("Oxy_density_neutron.txt");
    myproton.open("Oxy_density_proton.txt");
	
    int fileLine = 0;
    double radius, fileDensity;


    while(myneutron >> radius >> fileDensity)
	{
		density_neutron[fileLine] = fileDensity;
		fileLine++;
	}

    fileLine = 0;

    while(myproton >> radius >> fileDensity)
	{
		density_proton[fileLine] = fileDensity;
		fileLine++;
	}

   // for(int i=0; i<n_step_width_box; i++) cout << density_proton[i] << "\t" << density_neutron[i] << endl;
/*
    for(int i=0, i<n_step_width_box, i++){
    	density_neutron[i] << density_neutron << endl;
    }

    for(int i=0, i<n_step_width_box; i++){
    	density_proton[i] << density_proton << endl;
    }
*/    myneutron.close();
    myproton.close();



do{

	old_integral=integral;

	cout << "\nMy Integral: " << integral << "\tOld integral: " << old_integral << endl;

	cout << "\nPRESS ENTER TO CONTINUE!" << endl;
	cin.ignore(); 


// First: get the skyrme potential

    U_skyrme_p.clear();
	U_skyrme_n.clear();

	for(int i=0; i<n_step_width_box; i++){
		double r  = density_proton[i] + density_neutron[i];
		double rp = density_proton[i];
		double rn = density_neutron[i];
		U_skyrme_p.push_back(skyrme(r,rp,rp,rn));
		U_skyrme_n.push_back(skyrme(r,rn,rp,rn));
	}

//Routine to find eigevalues of neutron

    state appoggio;
    vector<double> wfWork;
	for(int node=0; node<3; node++){
		for(int L=0; L<8;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 0.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(0.05);
					wave_val.push_back(0.1);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, -1./2., U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
						if(wave_val[i-1]*wave_val[i]<0.)
						nodecount++;
					}
					if(nodecount>node){
						Eup = Etrial;
					}
					else if(nodecount<node||nodecount==node){
						Edown = Etrial;
					}
					wave_val.clear();
				}while(abs(Eup-Edown)>prec);
				if(Etrial<-0.1&&Etrial>-50){

					//Normalisation
					wfWork.push_back(0.05);
					wfWork.push_back(0.1);
					double norm = 0.;
					for(int i=1; i<n_step_width_box; i++) wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, -1./2.,  U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
					for(int i=1; i<n_step_width_box+1; i++) norm += h_width*pow(wfWork[i],2);
					wfWork.clear();
					//end normalisatio-1./2.

					appoggio.wavefn.push_back(0.05);
					appoggio.wavefn.push_back(0.1);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, -1./2., U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]/sqrt(norm)));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					neutronstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}


//This routine sort the vector of NEUTRON and print it

	sort(neutronstates.begin(), neutronstates.end(), compare());
	int n_neutron = 0, n_levels = 0;
	cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<neutronstates.size(); i++){
		if(n_neutron<126){
			cout<<neutronstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<neutronstates[i].eig<<endl;
			n_neutron += (2*neutronstates[i].j+1);
			n_levels ++;
		}
	}
	cout<<"Il numero di neutroni richiesto è:"<<n_neutron<<endl<<endl;


//Routine to find eigevalues of proton

    vector<state> protonstates;
	for(int node=0; node<3; node++){
		for(int L=0; L<8;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 10.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(0.05);
					wave_val.push_back(0.1);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
						if(wave_val[i-1]*wave_val[i]<0.)
						nodecount++;
					}
					if(nodecount>node){
						Eup = Etrial;
					}
					else if(nodecount<node||nodecount==node){
						Edown = Etrial;
					}
					wave_val.clear();
				}while(abs(Eup-Edown)>prec);
				if(Etrial<9.9&&Etrial>-50){

					//Normalisation
					wfWork.push_back(0.05);
					wfWork.push_back(0.1);
					double norm = 0.;
					for(int i=1; i<n_step_width_box; i++) wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, 1./2.,  U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1]));
					for(int i=1; i<n_step_width_box+1; i++) norm += h_width*pow(wfWork[i],2);
					wfWork.clear();
					//end normalisation

					appoggio.wavefn.push_back(0.05);
					appoggio.wavefn.push_back(0.1);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_p[i-1], U_skyrme_p[i], U_skyrme_p[i+1])/sqrt(norm));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					protonstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}

//This routine sort the vector of PROTON and print it

	sort(protonstates.begin(), protonstates.end(), compare());
	int n_proton = 0, p_levels = 0;
	cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<protonstates.size(); i++){
		if(n_proton<82){
			cout<<protonstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<protonstates[i].eig<<endl;
			n_proton += (2*protonstates[i].j+1);
			p_levels ++;
		}
	}
	cout<<"Il numero di protoni richiesto è:"<<n_proton<<endl;


//This routine calculate the new density for neutron

	for(int i=0; i<n_step_width_box+1; i++){
		density_neutron[i] = 0.;
	}
	for(int k=0; k<n_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_neutron[i] += (2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(neutronstates[k].wavefn[i],2);
	}

//This routine calculate the new density for proton

	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}
	for(int k=0; k<p_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	}


//////////////////TOTAL ENERGY///////////////////////////////////

//Calculate the kinetic density
	for(int i=0; i<n_step_width_box+1; i++){
		K_skyrme[i] = 0.;
	}
	for(int k=0; k<n_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			K_skyrme[i] += m_factor*(2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*(pow((neutronstates[k].wavefn[i]-neutronstates[k].wavefn[i-1])/h_width-neutronstates[k].wavefn[i]/(i*h_width+h_width),2)+
									neutronstates[k].l*(neutronstates[k].l+1)/pow(i*h_width+h_width,2)*pow(neutronstates[k].wavefn[i],2) );
	}

//Calculate the skyrme density

	for(int i=0; i<n_step_width_box+1; i++){
		E_skyrme[i] = 0.;
	}

	for(int i=2; i<n_step_width_box; i++){
		double r  = density_proton[i] + density_neutron[i];
		double rp = density_proton[i];
		double rn = density_neutron[i];
		E_skyrme[i] += 	1./2.*t0*((1+x0/2.)*pow(r,2)-(x0+1./2.)*(pow(rp,2)+pow(rn,2)))+ 
								1./12.*t3*pow(r,_a)*((1+x3/2.)*pow(r,2)-(x3+1./2.)*(pow(rp,2)+pow(rn,2)));
			cout<<"K_skyrme:"<<"\t"<<K_skyrme[i]<<"\t"<<"E_skyrme:"<<"\t"<<E_skyrme[i]<<endl;
	}



//Integral
	integral = 0.;
	for(int i=0; i<n_step_width_box+1; i++){
		integral += 4*M_PI*pow(i*h_width+h_width,2)*(K_skyrme[i]+E_skyrme[i])*h_width;
	}

///////////////////////////////////////////////////////////////

	cout << "\nMy Integral: " << integral << "\tOld integral: " << old_integral << endl;

	cout << "\nPRESS ENTER TO CONTINUE!" << endl;
	cin.ignore(); 

//cout<<abs(old_integral)<<endl;
}while(fabs(old_integral - integral) > prec);





//Routine to obtain proton density, and the sum of both
/*
	double * density_proton = new double [n_step_width_box+1];
	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}
	for(int k=0; k<p_levels-7; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	}

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

	*/





//USEFUL ROUTINES

//This is just for plotting all wavefunction
/*
	char orbital[26];
	ofstream myfile_wave;
	for(int k=0; k<p_levels-7; k++){
					sprintf(orbital,"wf_proton%1.4f.dat",protonstates[k].eig);
					myfile_wave.open(orbital);
		for(int i=2; i<n_step_width_box+1; i++){
					myfile_wave << i*h_width +h_width<< "   " << protonstates[k].wavefn[i]/(i*h_width +h_width) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
				}
				myfile_wave.close();
			}

*/

//This can be useful to plot a potential
/*
	ofstream myfile;
	myfile.open("Eigen.txt");

	for(int i=2; i<n_step_width_box; i++){
		wave_val.push_back(potential_woods(i*h_width +h_width)-centrifug_term(i*h_width +h_width, 0));
		myfile <<i*h_width +h_width<< "   " << wave_val[i]<< endl;
		myfile <<i*h_width +h_width<< "\t" << U_skyrme_p[i] << "\t" << U_skyrme_n[i] << endl;
	}

	myfile.close();
*/

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
