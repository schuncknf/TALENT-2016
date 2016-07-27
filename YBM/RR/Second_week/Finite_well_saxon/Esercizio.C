#include "Esercizio.h"

int main(){


//Finite well

	double Edown = 0.;
	double Eup = 100.;
	double n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	vector<double> eigenval;
	

//Routine to find eigevalues: this works finding the value of energy at which the eigenfunction cross the 0 
	
	double S = 1./2.;
	cout<<"Level"<<"    	"<<"L value"<<"    	"<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int node=0; node<4; node++){
		for(int L=0; L<8;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 0.;
				Edown = -100.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(0.);
					wave_val.push_back(0.05);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J));
					}
					nodecount = 0;
					for(int i=1; i<=n_step_width_box; i++){
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
			if(Etrial<-0.1&&Etrial>-80)
			cout<<node+1<<"	  	"<<L<<"    	 	"<<J*2<<"/2"<<"  	 		"<<Etrial<<endl;
			}
		}
		cout<<endl;
	}


// Write results to file with filename
/*
	ofstream myfile;
	myfile.open("Eigen.txt");
	int count =0;

	wave_val.push_back(0);
	wave_val.push_back(1E-2);
	for(int i=1; i<=n_step_width_box; i++){
		wave_val.push_back(numerov_algorithm_woods(eigenval[0], wave_val[i], wave_val[i-1], i*h_width+h_width,S,L));

		myfile << i*h_width +h_width<< "   " << wave_val[i] << endl;
	}

	myfile.close();
*/


//Plot the potential - Works
/*
	ofstream myfile;
	myfile.open("Eigen.txt");

	for(int i=0; i<n_step_width_box; i++){
		wave_val.push_back(potential_woods(i*h_width+2*h_width)+potential_spin_orbit(i*h_width+2*h_width,S, L)+potential_coulomb(i*h_width+2*h_width));
		myfile <<i*h_width +2*h_width<< "   " << wave_val[i] << endl;
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

	for(int i=0; i<n_step_width_box; i++){
		Proton_pot.push_back(potential_woods(i*h_width+2*h_width));
		Neutron_pot.push_back(potential_woods(i*h_width+2*h_width)+potential_coulomb(i*h_width+2*h_width));
		myfileN <<i*h_width +2*h_width<< "   " << Proton_pot[i] << endl;
		myfileP <<i*h_width +2*h_width<< "   " << Neutron_pot[i] << endl;
	}

	myfileN.close();
	myfileP.close();
*/

return 0;
}