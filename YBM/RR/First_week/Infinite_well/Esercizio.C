#include "Esercizio.h"

int main(){


//Infinite Well
//Discretize the energy to find eigenvalues
	double Edown = 0.;
	double Eup = 100.;
	double n_step_width = width/h_width;


	vector<double> wave_val;
	vector<double> eigenval;

//Find Eigenvalues

	for(int node=0; node<4; node++){
		int nodecount=0;
		double Etrial=0.;
		Edown = 0.;
		Eup = 100.;
		do{
			Etrial = (Eup+Edown)/2.;
			wave_val.push_back(0.);
			wave_val.push_back(0.05);
			for(int i=1; i<=n_step_width; i++){
				wave_val.push_back(numerov_algorithm(Etrial, wave_val[i], wave_val[i-1]));
			}

		nodecount = 0;
		for(int i=1; i<=n_step_width; i++){
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
	eigenval.push_back(Etrial);
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

	for(int i=1; i<=n_step_width; i++){
		wave_val0.push_back(numerov_algorithm(eigenval[0], wave_val0[i], wave_val0[i-1]));
		wave_val1.push_back(numerov_algorithm(eigenval[1], wave_val1[i], wave_val1[i-1]));
		wave_val2.push_back(numerov_algorithm(eigenval[2], wave_val2[i], wave_val2[i-1]));
		wave_val3.push_back(numerov_algorithm(eigenval[3], wave_val3[i], wave_val3[i-1]));
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