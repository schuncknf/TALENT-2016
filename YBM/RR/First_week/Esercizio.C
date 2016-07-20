#include "Esercizio.h"

int main(){


	double n_step_width = width/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;

	vector<double> wave_val;
	vector<double> eigenval;

//Discretising energy strategy to find eigenvalues

	double energy=0., wave_prev=0;

	for(int j=0; j<n_step_energy; j++){
		wave_val.push_back(0);
		wave_val.push_back(0.05);
		energy = Edown + j*h_energy;

		for(int i=0; i<n_step_width; i++){
			wave_val.push_back(numerov_algorithm(energy, wave_val[i+1], wave_val[i]));
		}

		if(wave_val[wave_val.size()-1]*wave_prev<0){
			cout<<energy-h_energy<<" "<<energy<<endl;
			eigenval.push_back(energy+h_energy/2.);
		}

		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();
	}

//Drawing the function

	

return 0;
}

		//cout<<energy<<"	"<<wave_val[wave_val.size()-1]<<endl;//<<"	"<<wave_prev<<endl;