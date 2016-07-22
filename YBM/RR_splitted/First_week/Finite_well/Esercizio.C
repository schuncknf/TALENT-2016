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
			wave_val.push_back(numerov_algorithm_finitewell(energy, wave_val[i+1], wave_val[i], i*h_width));
		}
		if(wave_val[wave_val.size()-1]*wave_prev<0.){
			cout<<energy-h_energy<<" "<<energy<<endl;
			eigenval.push_back(energy+h_energy/2.);
		}

		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();	
	}


		// Write results to file with filename including energy eigenvalue

	for(int j=0; j<n_step_energy; j++){
		wave_val.push_back(0.);
		wave_val.push_back(1E-5);
		energy = Edown + j*h_energy;

		for(int i=0; i<n_step_width_box; i++){
			wave_val.push_back(numerov_algorithm_finitewell(eigenvalue[j], wave_val[i+1], wave_val[i], i*h_width));
			opFile << wave_val.at(i) <<endl;
		}

		if(wave_val[wave_val.size()-1]*wave_prev<0.){
			cout<<energy-h_energy<<" "<<energy<<endl;
			eigenval.push_back(energy+h_energy/2.);
		}

		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();	
	}
				char filename[512], eigEng[512];
				sprintf(filename,"output_%1.4f.dat",eigenval[i]);
				sprintf(eigEng,"%1.4f",eigenEng);
				ofstream opFile;
				opFile.open(filename);
				for(int k=0; k<a/h; k++) opFile << h*k << "\t" << wf_val.at(k) << endl;
				opFile.close();


//Drawing the function
/*
	TApplication app("app", 0, 0);
	TGraph * wavefunction_finitewell = new TGraph();

	//This is to initialise the wavefunctions
	//The new vector that needs to subtract the function, probabily useless

		vector<double> wave_val2;
		wave_val.push_back(0.);
		wave_val.push_back(1E-5);	
		wave_val2.push_back(0.);
		wave_val2.push_back(1E-5);


	int count =0;


	for(int i=0; i<n_step_width_box; i++){		//Loop to draw the eigenfunction
		wave_val.push_back(numerov_algorithm_finitewell(eigenval[0], wave_val[i+1], wave_val[i], i*h_width));	//I use the energy I get after the first eigen
		if(i*h_width>20.&&wave_val[i]*wave_val[i-1]<0){
			count = i;
			break;
		}
		wavefunction_finitewell->SetPoint(i, i*h_width, wave_val[i]);
	}

	for(int i=count; i<n_step_width_box; i++){
		wave_val.push_back(0.);
		wavefunction_finitewell->SetPoint(i, i*h_width, wave_val[i]); //Fill the restant part of the wavefunctions with 0
	}

//Drawing potential
		
	for(int i=0; i<n_step_width_box; i++){
		wavefunction_finitewell->SetPoint(i, i*h_width, potential_woods(i*h_width-12)); //Fill the restant part of the wavefunctions with 0
	}
	

		TCanvas myCanvas("tela","tela");
		wavefunction_finitewell->SetLineColor(2);
	   	wavefunction_finitewell->GetXaxis()->SetTitle("l [fm]");
   		wavefunction_finitewell->GetYaxis()->SetTitle("#Psi");
   		//wavefunction_finitewell->GetYaxis()->SetRangeUser(-40.,40.); 
   	 	wavefunction_finitewell->Draw("AC");

    app.Run();
	delete wavefunction_finitewell;
*/
return 0;
}