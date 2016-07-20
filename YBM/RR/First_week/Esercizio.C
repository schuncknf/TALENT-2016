#include "Esercizio.h"

int main(){


	double n_step_width = width/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;


	vector<double> wave_val;
	wave_val.push_back(0);
	wave_val.push_back(0.05);


//Discretising energy strategy

	double energy=0., wave_last, wave_prev=0;
	//for(int j=0; j<n_step_energy; j++){
		//energy = Edown + j*h_energy;
		energy = 0.668;
		for(int i=0; i<n_step_width; i++){
			wave_val.push_back(numerov_algorithm(energy, wave_val[i+1], wave_val[i]));
			//cout<<wave_val[i]<<endl;
		}
		cout<<energy<<"	"<<wave_val[wave_val.size()-1]<<endl;//<<"	"<<wave_prev<<endl;
		if(wave_val[wave_val.size()-1]*wave_prev<0){
			//cout<<energy-h_energy<<" "<<energy<<endl;
		}
		wave_prev = wave_val[wave_val.size()-1];
		//cout<<wave_prev<<endl;
		wave_val.clear();
		//cout<<energy<<endl;
	//}

return 0;
}
