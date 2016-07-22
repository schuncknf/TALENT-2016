#include "Esercizio.h"

int main(){


//Finite well
//Discretize the energy to find eigenvalues

	double n_step_width_box = width_box/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;
	double energy=0., wave_prev=0;

	vector<double> wave_val;
	vector<double> eigenval;
	

	//Routine to find eigevalues: this works finding the value of energy at which the eigenfunction cross the 0 

	for(int j=0; j<n_step_energy; j++){
		wave_val.push_back(0.);
		wave_val.push_back(1E-5);
		energy = Edown + j*h_energy;

		for(int i=0; i<n_step_width_box; i++){
			wave_val.push_back(numerov_algorithm_woods(energy, wave_val[i+1], wave_val[i], i*h_width));
		}
		if(wave_val[wave_val.size()-1]*wave_prev<0.){
			cout<<energy-h_energy<<" "<<energy<<endl;
			eigenval.push_back(energy+h_energy/2.);
		}

		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();	
	}


// Write results to file with filename

	ofstream myfile;
	myfile.open("Eigen.txt");
	int count =0;

	wave_val.push_back(0);
	wave_val.push_back(1E-5);
	for(int i=0; i<n_step_width_box+1; i++){
		wave_val.push_back(numerov_algorithm_woods(eigenval[3], wave_val[i+1], wave_val[i], i*h_width));
		if(i*h_width>24.&&wave_val[i]*wave_val[i-1]<0){
			count = i;
			break;
		}
		myfile << i*h_width << "   " << wave_val[i] << endl;
	}
	for(int i=count; i<n_step_width_box; i++){
		wave_val.push_back(0.);
		myfile << i*h_width << "   " << wave_val[i] << endl;
	}


	myfile.close();

return 0;
}