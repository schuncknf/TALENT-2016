#include "Esercizio.h"

int main(){


//Finite well

	double Edown = 0.;
	double Eup = 100.;
	double n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	vector<double> eigenval;
	

//Routine to find eigevalues: this works finding the value of energy at which the eigenfunction cross the 0 

	for(int node=0; node<4; node++){
		int nodecount=0;
		double Etrial=0.;
		Edown = 0.;
		Eup = 100.;
		do{
			Etrial = (Eup+Edown)/2.;
			wave_val.push_back(0.);
			wave_val.push_back(0.05);
			for(int i=1; i<=n_step_width_box; i++){
				wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width));
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
	eigenval.push_back(Etrial);
}


// Write results to file with filename

	ofstream myfile;
	myfile.open("Eigen.txt");
	int count =0;

	wave_val.push_back(0);
	wave_val.push_back(1E-5);
	for(int i=0; i<=n_step_width_box; i++){
		wave_val.push_back(numerov_algorithm_woods(eigenval[n], wave_val[i+1], wave_val[i], i*h_width));
		if(i*h_width>24.&&wave_val[i]*wave_val[i-1]<0){
			count = i;
			break;
		}
		myfile << i*h_width << "   " << wave_val[i] << endl;
	}
	for(int i=count; i<=n_step_width_box; i++){
		wave_val.push_back(0.);
		myfile << i*h_width << "   " << wave_val[i] << endl;
	}

	cout<<"The eigenvalues are:"<<endl;
	for(int i=0; i<eigenval.size(); i++){
		cout<<eigenval[i]<<endl;
	}	


	myfile.close();

return 0;
}