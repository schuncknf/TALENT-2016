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

	cout<<"Which eigenfunction do you want?"<<endl;
	int a;	cin>>a;
	TApplication app("app", 0, 0);
	TGraph * wavefunction = new TGraph();
		wave_val.push_back(0);
		wave_val.push_back(0.05);
	for(int i=0; i<n_step_width; i++){
		wave_val.push_back(numerov_algorithm(eigenval[a], wave_val[i+1], wave_val[i]));
		wavefunction->SetPoint(i, i*h_width, wave_val[i]);
	}
	TCanvas myCanvas("tela","tela");
	wavefunction->SetLineColor(2);
	wavefunction->GetXaxis()->SetRangeUser(0.,35.); 
    wavefunction->GetXaxis()->SetTitle("l [fm]");
    wavefunction->GetYaxis()->SetTitle("#Psi");
    wavefunction->Draw("AC");

    app.Run();
	delete wavefunction;


return 0;
}