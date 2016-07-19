#include "Esercizio.h"

int main(){

	double a = 35;
	double h_width = 0.01;
	double n_step_width = a/h_width;

	double h_energy = 0.0005;
	double n_step_energy = (Eup-Edown)/h_energy;


	vector<double> wave_val;
	wave_val.push_back(0);
	wave_val.push_back(0.5);


//Discretising energy strategy

	double energy, wave_last, wave_prev=0;
	for(int j=0; j<n_step_energy; j++){
		energy = Edown + j*h_energy;

		for(int i=0; i<n_step_width; i++){
			wave_val.push_back(numerov_algorithm(energy, h_width, wave_val[i+1], wave_val[i]));
			//cout<<wave_val[i]<<endl;
		}
		if(wave_val[wave_val.size()-1]*wave_prev<0){
			cout<<energy-h_energy<<" "<<energy<<endl;
		}
		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();
	}

return 0;
}
