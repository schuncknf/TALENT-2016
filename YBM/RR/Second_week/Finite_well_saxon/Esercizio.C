#include "Esercizio.h"

int main(){


//Finite well

	double Edown = 0.;
	double Eup = -50.;
	int n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	vector<double> eigenval;
	double S = 1./2.;

//Routine to find eigevalues and write result on terminal


	int n_neutron = 0;
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
					wave_val.push_back(0.05);
					wave_val.push_back(0.1);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,-1./2.));
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
					n_neutron += (2*J+1);
					cout<<node+1<<"	  	"<<L<<"    	 	"<<J*2<<"/2"<<"  	 		"<<setprecision(12)<<Etrial<<endl;
				}	
			}
		}
		cout<<endl;
	}
	cout<<"The number of neutron is:"<<n_neutron<<endl;



//Routine to obtain neutron density and proton density, and the sum of both
/*
	ofstream myfile_density, myfile_wave;
	char orbital[26];
	int n_neutron = 0;

	double * density = new double [n_step_width_box+1];
	for(int i=1; i<n_step_width_box+1; i++){
		density[i] = 0.;
	}
	vector<double> wave_eigen;
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
					for(int i=1; i<n_step_width_box; i++){ //the vector has n_step_width_box + 1 boxes
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,-1./2.));
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
				if(Etrial<-4.5){
					wave_eigen.push_back(0.05);
					wave_eigen.push_back(0.1);
					sprintf(orbital,"wf_neutron%1.4f.dat",Etrial);
					myfile_wave.open(orbital);
					//cout<<Etrial<<endl;
					double norm = normalise(Etrial, n_step_width_box, S, L, J,-1./2.);
					for(int i=1; i<n_step_width_box+1; i++){
						wave_eigen.push_back(numerov_algorithm_woods(Etrial, wave_eigen[i], wave_eigen[i-1], i*h_width+h_width,S,L,J,-1./2.));
						density[i] += (2.*J+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(wave_eigen[i]/sqrt(norm),2);
						myfile_wave << i*h_width +h_width<< "   " << wave_eigen[i]/(sqrt(norm)*(i*h_width +h_width)) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
					}
					wave_eigen.clear();
					myfile_wave.close();
				}
			}
		}	
	}
	myfile_density.open("density_neutron.txt");
	for(int i=1; i<n_step_width_box+1; i++){
		myfile_density << i*h_width +h_width<< "   " << density[i] << endl;
	}
	myfile_density.close();



	ofstream myfile_density_proton, myfile_wave_proton;
	char orbital_proton[26];

	double * density_proton = new double [n_step_width_box+1];
	for(int i=1; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}
	vector<double> wave_eigen_proton;
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
					for(int i=1; i<n_step_width_box; i++){ //the vector has n_step_width_box + 1 boxes
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,1./2.));
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
				if(Etrial<8.4&&Etrial>-50){
					wave_eigen_proton.push_back(0.05);
					wave_eigen_proton.push_back(0.1);
					sprintf(orbital_proton,"wf_proton%1.4f.dat",Etrial);
					myfile_wave_proton.open(orbital_proton);
					//cout<<Etrial<<"	"<<endl;
					double norm = normalise(Etrial, n_step_width_box, S, L, J,1./2.);
					for(int i=1; i<n_step_width_box+1; i++){
						wave_eigen_proton.push_back(numerov_algorithm_woods(Etrial, wave_eigen_proton[i], wave_eigen_proton[i-1], i*h_width+h_width,S,L,J,1./2.));
						density_proton[i] += (2.*J+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(wave_eigen_proton[i]/sqrt(norm),2);
						myfile_wave_proton << i*h_width +h_width<< "   " << wave_eigen_proton[i]/(sqrt(norm)*(i*h_width +h_width)) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
					}
					wave_eigen_proton.clear();
					myfile_wave_proton.close();
				}
			}
		}	
	}

	myfile_density_proton.open("density_proton.txt");
	for(int i=1; i<n_step_width_box+1; i++){
		myfile_density_proton << i*h_width +h_width<< "   " << density_proton[i] << endl;
	}

	myfile_density_proton.close();

	ofstream myfile_density_all;
	myfile_density_all.open("density_all.txt");
	for(int i=1; i<n_step_width_box+1; i++){
		myfile_density_all << i*h_width +h_width<< "   " << density_proton[i]+density[i] << endl;
	}
	myfile_density_all.close();
*/

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


//Plot a potential
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