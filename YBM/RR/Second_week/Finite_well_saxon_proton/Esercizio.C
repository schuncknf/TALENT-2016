#include "Esercizio.h"

int main(){


//Finite well

	double Edown = 0.;
	double Eup = -50.;
	int n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	vector<double> eigenval;
	

//Routine to find eigevalues: this works finding the value of energy at which the eigenfunction cross the 0 
/*	
	double S = 1./2.;
	cout<<"Level"<<"    	"<<"L value"<<"    	"<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
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
					wave_val.push_back(0.1);
					wave_val.push_back(0.15);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+2.*h_width,S,L,J));
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
			if(Etrial<-4.5&&Etrial>-50)
			cout<<node+1<<"	  	"<<L<<"    	 	"<<J*2<<"/2"<<"  	 		"<<setprecision(12)<<Etrial<<endl;
			}
		}
		cout<<endl;
	}
*/


//Density
	ofstream myfile_density, myfile_wave;
	char orbital[26];
	//myfile_density.open("density.txt");

	double * density = new double [n_step_width_box+1];
	for(int i=1; i<n_step_width_box+1; i++){
		density[i] = 0.;
	}
	vector<double> wave_eigen;
	double S = 1./2.;
	for(int node=0; node<3; node++){
	//int node =0;
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
					wave_val.push_back(0.1);
					wave_val.push_back(0.15);
					for(int i=1; i<n_step_width_box; i++){ //the vector has n_step_width_box + 1 boxes
						wave_val.push_back(numerov_algorithm_woods_proton(Etrial, wave_val[i], wave_val[i-1], i*h_width+2.*h_width,S,L,J));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+1; i++){
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
				if(Etrial<-4.5&&Etrial>-50){
					wave_eigen.push_back(0.1);
					wave_eigen.push_back(0.15);
					sprintf(orbital,"wf_%1.4f_.dat",Etrial);
					myfile_wave.open(orbital);
					cout<<Etrial<<"	";
					double norm = normalise(Etrial, n_step_width_box, S, L, J);
					cout <<norm<<endl;
					for(int i=1; i<n_step_width_box+1; i++){
						wave_eigen.push_back(numerov_algorithm_woods_proton(Etrial, wave_eigen[i], wave_eigen[i-1], i*h_width+2.*h_width,S,L,J));
						density[i] += (2.*J+1.)/(4.*M_PI*pow(i*h_width+2.*h_width,2))*pow(wave_eigen[i]/sqrt(norm),2);
						myfile_wave << i*h_width +2.*h_width<< "   " << wave_eigen[i]/(sqrt(norm)*(i*h_width +2.*h_width)) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
					}
					wave_eigen.clear();
					myfile_wave.close();
				}
			}
		}	
	}

	myfile_density.open("Density.txt");
	for(int i=1; i<n_step_width_box+1; i++){
		myfile_density << i*h_width +2.*h_width<< "   " << density[i] << endl;
	}


	myfile_density.close();

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

//a=system('a=`tempfile`;cat *.dat > $a;echo "$a"') 
//Plot the potential - Works
/*
	ofstream myfile;
	myfile.open("Eigen.txt");

	for(int i=0; i<n_step_width_box; i++){
		wave_val.push_back(potential_woods(i*h_width +h_width)-centrifug_term(i*h_width +h_width, 0));
		myfile <<i*h_width +h_width<< "   " << wave_val[i] << endl;
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