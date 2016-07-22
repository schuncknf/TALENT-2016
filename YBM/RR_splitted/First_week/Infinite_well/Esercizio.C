#include "Esercizio.h"

int main(){


//Infinite Well
//Discretize the energy to find eigenvalues

	double n_step_width = width/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;

	vector<double> wave_val;
	vector<double> eigenval;

	double energy=0., wave_prev=0;

//Find Eigenvalues

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

//Write eigenfunctions data on files

	ofstream myfile0,myfile1,myfile2,myfile3;
	myfile0.open("Eigen_0.txt");
	myfile1.open("Eigen_1.txt");
	myfile2.open("Eigen_2.txt");
	myfile3.open("Eigen_3.txt");

	vector<double> wave_val0;
	vector<double> wave_val1;
	vector<double> wave_val2;
	vector<double> wave_val3;

	wave_val0.push_back(0);
	wave_val0.push_back(0.05);
	wave_val1.push_back(0);
	wave_val1.push_back(0.05);
	wave_val2.push_back(0);
	wave_val2.push_back(0.05);
	wave_val3.push_back(0);
	wave_val3.push_back(0.05);

	for(int i=0; i<n_step_width+1; i++){
		wave_val0.push_back(numerov_algorithm(eigenval[0], wave_val0[i+1], wave_val0[i]));
		wave_val1.push_back(numerov_algorithm(eigenval[1], wave_val1[i+1], wave_val1[i]));
		wave_val2.push_back(numerov_algorithm(eigenval[2], wave_val2[i+1], wave_val2[i]));
		wave_val3.push_back(numerov_algorithm(eigenval[3], wave_val3[i+1], wave_val3[i]));
		myfile0 << i*h_width << "   " << wave_val0[i] << endl;
		myfile1 << i*h_width << "   " << wave_val1[i] << endl;
		myfile2 << i*h_width << "   " << wave_val2[i] << endl;
		myfile3 << i*h_width << "   " << wave_val3[i] << endl;
	}

	myfile0.close();
	myfile1.close();
	myfile2.close();
	myfile3.close();

return 0;
}