#include "Esercizio.h"

int main(){


//Infinite Well
//Discretize the energy to find eigenvalues
/*
	double n_step_width = width/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;

	vector<double> wave_val;
	vector<double> eigenval;

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
	vector<double> wave_val0;
	vector<double> wave_val1;
	vector<double> wave_val2;
	vector<double> wave_val3;
	vector<double> wave_val4;
	vector<double> wave_val5;

	TApplication app("app", 0, 0);
	TGraph * wavefunction_0 = new TGraph();
	TGraph * wavefunction_1 = new TGraph();
	TGraph * wavefunction_2 = new TGraph();
	TGraph * wavefunction_3 = new TGraph();
	TGraph * wavefunction_4 = new TGraph();
	TGraph * wavefunction_5 = new TGraph();

	wave_val0.push_back(0);
	wave_val0.push_back(0.05);
	wave_val1.push_back(0);
	wave_val1.push_back(0.05);
	wave_val2.push_back(0);
	wave_val2.push_back(0.05);
	wave_val3.push_back(0);
	wave_val3.push_back(0.05);
	wave_val4.push_back(0);
	wave_val4.push_back(0.05);
	wave_val5.push_back(0);
	wave_val5.push_back(0.05);
	for(int i=0; i<n_step_width; i++){
		wave_val0.push_back(numerov_algorithm(eigenval[0], wave_val0[i+1], wave_val0[i]));
		wave_val1.push_back(numerov_algorithm(eigenval[1], wave_val1[i+1], wave_val1[i]));
		wave_val2.push_back(numerov_algorithm(eigenval[2], wave_val2[i+1], wave_val2[i]));
		wave_val3.push_back(numerov_algorithm(eigenval[3], wave_val3[i+1], wave_val3[i]));
		wave_val4.push_back(numerov_algorithm(eigenval[4], wave_val4[i+1], wave_val4[i]));
		wave_val5.push_back(numerov_algorithm(eigenval[5], wave_val5[i+1], wave_val5[i]));

		wavefunction_0->SetPoint(i, i*h_width, wave_val0[i]);
		wavefunction_1->SetPoint(i, i*h_width, wave_val1[i]);
		wavefunction_2->SetPoint(i, i*h_width, wave_val2[i]);
		wavefunction_3->SetPoint(i, i*h_width, wave_val3[i]);
		wavefunction_4->SetPoint(i, i*h_width, wave_val4[i]);
		wavefunction_5->SetPoint(i, i*h_width, wave_val5[i]);
	}

	TCanvas myCanvas("tela","tela");
    wavefunction_0->SetLineColor(2);
    wavefunction_1->SetLineColor(4);
    wavefunction_2->SetLineColor(9);
    wavefunction_3->SetLineColor(8);
    wavefunction_4->SetLineColor(6);
    wavefunction_5->SetLineColor(1);
  	wavefunction_0->GetXaxis()->SetRangeUser(0.,35.); 
  	wavefunction_0->GetYaxis()->SetRangeUser(-30.,60.); 
    wavefunction_0->GetXaxis()->SetTitle("l [fm]");
    wavefunction_0->GetYaxis()->SetTitle("#Psi");

    wavefunction_0->Draw("AC");
    wavefunction_1->Draw("SAME");
    wavefunction_2->Draw("SAME");
    wavefunction_3->Draw("SAME");
    wavefunction_4->Draw("SAME");
    wavefunction_5->Draw("SAME");
    
    app.Run();
	delete wavefunction_0, wavefunction_1, wavefunction_2, wavefunction_3, wavefunction_4, wavefunction_5;
*/

//Finite well
//Discretize the energy to find eigenvalues

	double n_step_width_box = width_box/h_width;
	double n_step_energy = (Eup-Edown)/h_energy;
	double energy=0., wave_prev=0;

	vector<double> wave_val;
	vector<double> eigenval;

	for(int j=0; j<n_step_energy; j++){
		wave_val.push_back(0.);
		wave_val.push_back(1E-5);
		energy = Edown + j*h_energy;

		for(int i=0; i<n_step_width_box; i++){
			wave_val.push_back(numerov_algorithm_finitewell(energy, wave_val[i+1], wave_val[i], i*h_width));
		}
		if(wave_val[wave_val.size()-1]*wave_prev<0){
			cout<<energy-h_energy<<" "<<energy<<endl;
			//eigenval.push_back(energy+h_energy/2.);
		}

		wave_prev = wave_val[wave_val.size()-1];
		wave_val.clear();	
	}



	//Drawing the function

	TApplication app("app", 0, 0);
	TGraph * wavefunction_finitewell = new TGraph();

		wave_val.push_back(0.);
		wave_val.push_back(1E-5);
	for(int i=0; i<n_step_width_box; i++){
		wave_val.push_back(numerov_algorithm_finitewell(7.2525, wave_val[i+1], wave_val[i], i*h_width));
		wavefunction_finitewell->SetPoint(i, i*h_width, wave_val[i]);
	}
		TCanvas myCanvas("tela","tela");
		wavefunction_finitewell->SetLineColor(2);
	   	wavefunction_finitewell->GetXaxis()->SetTitle("l [fm]");
   		wavefunction_finitewell->GetYaxis()->SetTitle("#Psi");
   	 	wavefunction_finitewell->Draw("AC");

    app.Run();
	delete wavefunction_finitewell;

return 0;
}