#include "Esercizio.h"

int main(){


//////////////////////////////////This first part needs to fill the vector with densities with the woods saxon////////////////////////////////////////
//Some costants and variables

	double Edown = 0.;
	double Eup = -50.;
	int n_step_width_box = width_box/h_width;

	vector<double> wave_val;
	vector<double> wfWork;
    vector<state> neutronstates, protonstates;
    vector<double> U_skyrme_p, U_skyrme_n, U_cou_p;
   	double integral_kin=0., integral_en=0., old_integral=1E5, S = 1./2.;
   	int n_levels, p_levels;

	double * density_neutron = new double [n_step_width_box];
	double * density_proton = new double [n_step_width_box];
	double * K_skyrme_n = new double [n_step_width_box];
	double * K_skyrme_p = new double [n_step_width_box];
	double * E_skyrme = new double [n_step_width_box];


    int fileLine = 0;
    double radius, fileDensity;
    double r1, rp1, rn1, r2, rp2, rn2, r3, rp3, rn3;
    double r, rp, rn;
   	double integ1 = 0., integ2 = 0., integ3 = 0., exchange = 0.;

//I get the first density from files

    state appoggio;
	for(int node=0; node<max_n; node++){
		for(int L=0; L<max_l;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 0.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(first_point);
					wave_val.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,-1./2.));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
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
				if(Etrial<-0.1&&Etrial>-50){
					double norm = normalise(Etrial, n_step_width_box, S, L, J,-1./2.);
					appoggio.wavefn.push_back(first_point);
					appoggio.wavefn.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,-1./2.)/sqrt(norm));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					neutronstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}


//This routine sort the vector of NEUTRON and print it

	sort(neutronstates.begin(), neutronstates.end(), compare());
	int n_neutron = 0; n_levels = 0;
	for(int i=0; i<neutronstates.size(); i++){
		if(n_neutron<N){
			n_neutron += (2*neutronstates[i].j+1);
			n_levels ++;
		}
	}


//Routine to find eigevalues of proton

    protonstates.clear();
	for(int node=0; node<max_n; node++){
		for(int L=0; L<max_l;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 10.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(first_point);
					wave_val.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,1./2.));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
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
				if(Etrial<9.9&&Etrial>-50){
					double norm = normalise(Etrial, n_step_width_box, S, L, J,1./2.);
					appoggio.wavefn.push_back(first_point);
					appoggio.wavefn.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_woods(Etrial, wave_val[i], wave_val[i-1], i*h_width+h_width,S,L,J,1./2.)/sqrt(norm));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					protonstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}

//This routine sort the vector of PROTON and print it

	sort(protonstates.begin(), protonstates.end(), compare());
	int n_proton = 0; p_levels = 0;
	for(int i=0; i<protonstates.size(); i++){
		if(n_proton<Z){
			n_proton += (2*protonstates[i].j+1);
			p_levels ++;
		}
	}

//Routine to obtain neutron density and proton density, and the sum of both

	for(int i=0; i<n_step_width_box+1; i++){
		density_neutron[i] = 0.;
	}
	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}

	for(int k=0; k<n_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_neutron[i] += (2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(neutronstates[k].wavefn[i],2);
	}

	for(int k=0; k<p_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	}



/////////////////////////////////////////////Here starts the HF loop////////////////////////////////////////////

int iteration =0;

//do{
//for(int u=0; u<2; u++){
 	old_integral = integral_en+integral_kin;

// First: get the skyrme potential

    U_skyrme_p.clear();
	U_skyrme_n.clear();
	U_cou_p.clear();

	for(int i=2; i<n_step_width_box-1; i++){

		//Skyrme potential
		U_skyrme_p.push_back(skyrme(density_proton, density_neutron, density_proton, i, i*h_width+h_width, K_skyrme_n[i] + K_skyrme_p[i], K_skyrme_p[i]));
		U_skyrme_n.push_back(skyrme(density_neutron, density_neutron, density_proton, i, i*h_width+h_width, K_skyrme_n[i] + K_skyrme_p[i], K_skyrme_n[i]));

		//Coulomb potential
		/*integ1 = 0., integ2 = 0., integ3 = 0., exchange = 0.;
		for(int count=0;count<i+1;count++){
				integ1 += density_proton[count]*pow(count*h_width+h_width,2);
				integ2 += density_proton[count]*(count*h_width+h_width);
		}
		for(int count=0; count<n_step_width_box+1;count++){
				integ3 += density_proton[count]*(count*h_width+h_width);
		}
		exchange = pow(density_proton[i],1./3.);
		U_cou_p.push_back( h_width* (4.*M_PI*e* (integ1/(i*h_width+h_width) - integ2 + integ3) - exchange*e*pow(3./M_PI,1./3.)) );
	*/}

	neutronstates.clear();
	protonstates.clear();

//Routine to find eigevalues of neutron

    appoggio.wavefn.clear();
    wfWork.clear();
	for(int node=0; node<max_n; node++){
		for(int L=0; L<max_l;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 0.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(first_point);
					wave_val.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, -1./2., U_skyrme_n[i-1], U_skyrme_n[i], U_skyrme_n[i+1]));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
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
				if(Etrial<-0.1&&Etrial>-50){
					//Normalisation
					wfWork.push_back(first_point);
					wfWork.push_back(second_point);
					double norm = 0.;
					for(int i=1; i<n_step_width_box+1; i++) {
						wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, -1./2.,  U_skyrme_n[i-1], U_skyrme_n[i], U_skyrme_n[i+1]));
					 	norm += h_width*pow(wfWork[i],2);
					}
					wfWork.clear();
					//end normalisation

					appoggio.wavefn.push_back(first_point);
					appoggio.wavefn.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, -1./2., U_skyrme_n[i-1], U_skyrme_n[i], U_skyrme_n[i+1])/sqrt(norm));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					neutronstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}


//This routine sort the vector of NEUTRON and print it

	sort(neutronstates.begin(), neutronstates.end(), compare());
	n_neutron = 0; n_levels = 0;
	//cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<neutronstates.size(); i++){
		if(n_neutron<N){
			//cout<<neutronstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<neutronstates[i].eig<<endl;
			n_neutron += (2*neutronstates[i].j+1);
			n_levels ++;
		}
	}
	//cout<<"Requested number of neutrons:"<<n_neutron<<endl;


//Routine to find eigevalues of proton

	for(int node=0; node<max_n; node++){
		for(int L=0; L<max_l;L++){
			for(int h=1; h<3; h++){
				if(L==0) h=2;
				double J = abs(S+pow(-1,h)*L);
				int nodecount=0;
				double Etrial=0.;
				Eup = 0.;
				Edown = -50.;
				do{
					Etrial = (Eup+Edown)/2.;
					wave_val.push_back(first_point);
					wave_val.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						wave_val.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_n[i-1] /*+ U_cou_p[i-1] */, U_skyrme_n[i]/*+ U_cou_p[i] */, U_skyrme_n[i+1]/*+ U_cou_p[i+1] */));
					}
					nodecount = 0;
					for(int i=1; i<n_step_width_box+2; i++){
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
				if(Etrial<9.9&&Etrial>-50){

					
					//Normalisation
					wfWork.push_back(first_point);
					wfWork.push_back(second_point);
					double norm = 0.;
					for(int i=1; i<n_step_width_box; i++){
						wfWork.push_back(numerov_algorithm_HF(Etrial, wfWork[i], wfWork[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_n[i-1] /*+ U_cou_p[i-1] */, U_skyrme_n[i]/*+ U_cou_p[i] */, U_skyrme_n[i+1]/*+ U_cou_p[i+1] */));
						norm += h_width*pow(wfWork[i],2);
					}
					wfWork.clear();
					//end normalisation
					

					appoggio.wavefn.push_back(first_point);
					appoggio.wavefn.push_back(second_point);
					for(int i=1; i<n_step_width_box; i++){
						appoggio.wavefn.push_back(numerov_algorithm_HF(Etrial, wave_val[i], wave_val[i-1], i*h_width + h_width, S, L, J, 1./2., U_skyrme_n[i-1] /*+ U_cou_p[i-1] */, U_skyrme_n[i]/*+ U_cou_p[i] */, U_skyrme_n[i+1]/*+ U_cou_p[i+1] */)/sqrt(norm));
					}
					appoggio.eig = Etrial;
					appoggio.l = L;
					appoggio.j = J;
					protonstates.push_back(appoggio);
					appoggio.wavefn.clear();
				}	
			}
		}
	}

//This routine sort the vector of PROTON and print it

	sort(protonstates.begin(), protonstates.end(), compare());
	n_proton = 0; p_levels = 0;
	cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<protonstates.size(); i++){
		if(n_proton<Z){
			cout<<protonstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<protonstates[i].eig<<endl;
			n_proton += (2*protonstates[i].j+1);
			p_levels ++;
		}
	}
	//cout<<"Requested number of protons:"<<n_proton<<endl<<endl;


//This routine calculate the new density for neutron

	for(int i=0; i<n_step_width_box+1; i++){
		density_neutron[i] = 0.;
	}
	for(int k=0; k<n_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_neutron[i] += (2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(neutronstates[k].wavefn[i],2);
	}

//This routine calculate the new density for proton

	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}
	for(int k=0; k<p_levels; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	}


//////////////////TOTAL ENERGY///////////////////////////////////

//Calculate the kinetic density
	for(int i=0; i<n_step_width_box+1; i++){
		K_skyrme_n[i] = 0.;
		K_skyrme_p[i] = 0.;
	}
	for(int k=0; k<n_levels; k++){
		for(int i=2; i<n_step_width_box; i++){
			K_skyrme_n[i] +=	m_factor*(1-1./A)*(2.*neutronstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*(pow((neutronstates[k].wavefn[i]-neutronstates[k].wavefn[i-1])/h_width-neutronstates[k].wavefn[i]/(i*h_width+h_width),2)+
									neutronstates[k].l*(neutronstates[k].l+1)/pow(i*h_width+h_width,2)*pow(neutronstates[k].wavefn[i],2));
		}
	}

	for(int k=0; k<p_levels; k++){
		for(int i=2; i<n_step_width_box; i++){
			K_skyrme_p[i] +=	m_factor*(1-1./A)*(2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*(pow((protonstates[k].wavefn[i]-protonstates[k].wavefn[i-1])/h_width-protonstates[k].wavefn[i]/(i*h_width+h_width),2)+
									protonstates[k].l*(protonstates[k].l+1)/pow(i*h_width+h_width,2)*pow(protonstates[k].wavefn[i],2) );
		}
	}

//Calculate the skyrme density

	for(int i=0; i<n_step_width_box+1; i++){
		E_skyrme[i] = 0.;
	}

	for(int i=2; i<n_step_width_box; i++){
		r  = density_proton[i] + density_neutron[i];
		rp = density_proton[i];
		rn = density_neutron[i];
		E_skyrme[i] += 	1./2.*t0*((1+x0/2.)*pow(r,2)-(x0+1./2.)*(pow(rp,2)+pow(rn,2)))+ 
						1./12.*t3*pow(r,_a)*((1+x3/2.)*pow(r,2)-(x3+1./2.)*(pow(rp,2)+pow(rn,2)));
	}



//Integral
	integral_kin = 0., integral_en= 0.;
	for(int i=2; i<n_step_width_box+1; i++){
		integral_kin += 4*M_PI*pow(i*h_width+h_width,2)*(K_skyrme_n[i]+K_skyrme_p[i])*h_width;
		integral_en += 4*M_PI*pow(i*h_width+h_width,2)*(E_skyrme[i])*h_width;
	}

cout<<"Kinetic energy:"<<integral_kin<<"\t"<<"Skyrme energy:"<<integral_en<<endl;
///////////////////////////////////////////////////////////////

iteration++;
//}
//}while(fabs(old_integral - integral_kin -integral_en) > HF_prec);

	cout<<"Neutrons:"<<endl;
	cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<p_levels; i++){
		cout<<protonstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<protonstates[i].eig<<endl;
	}
	cout<<"Protons:"<<endl;
	cout<<"J value"<<"  	  	"<<"Eigenvalue"<<endl;
	for(int i=0; i<n_levels; i++){
		cout<<neutronstates[i].j*2<<"/2"<<"  	 		"<<setprecision(12)<<neutronstates[i].eig<<endl;
	}

cout<<"Number of iteration:"<<iteration<<endl;
cout<<"Total energy:"<<integral_kin+integral_en<<endl;
cout<<"Energy difference:"<<abs(old_integral - integral_kin - integral_en)<<endl;



//Routine to obtain proton density, and the sum of both
/*
	double * density_proton = new double [n_step_width_box+1];
	for(int i=0; i<n_step_width_box+1; i++){
		density_proton[i] = 0.;
	}
	for(int k=0; k<p_levels-7; k++){
		for(int i=2; i<n_step_width_box; i++)
			density_proton[i] += (2.*protonstates[k].j+1.)/(4.*M_PI*pow(i*h_width+h_width,2))*pow(protonstates[k].wavefn[i],2);
	}

	ofstream myfile_density_neutron, myfile_density_proton, myfile_density_np;
	myfile_density_neutron.open("density_neutron.txt");
	myfile_density_proton.open("density_proton.txt");
	myfile_density_np.open("density_np.txt");

	for(int i=2; i<n_step_width_box+1; i++){
		myfile_density_neutron << i*h_width +h_width<< "   " << density_neutron[i] << endl;
		myfile_density_proton << i*h_width +h_width<< "   " << density_proton[i] << endl;
		myfile_density_np << i*h_width +h_width << "   " << density_neutron[i] + density_proton[i] << endl;
	}
	myfile_density_neutron.close();
	myfile_density_proton.close();
	myfile_density_np.close();

	*/


//USEFUL ROUTINES

//This is just for plotting all wavefunction
/*
cout<<n_levels<<endl;
	char orbital[26];
	ofstream myfile_wave;
	for(int k=0; k<n_levels; k++){
					sprintf(orbital,"wf_neutron%1.4f.dat",neutronstates[k].eig);
					myfile_wave.open(orbital);
		for(int i=2; i<n_step_width_box+1; i++){
					myfile_wave << i*h_width +h_width<< "   " << neutronstates[k].wavefn[i]/(i*h_width +h_width) << endl;	//Crea un numero di file pari agli autovalori e guarda nelle f d'onda dove ci sono i problemi ;-)
				}
				myfile_wave.close();
			}*/



//This can be useful to plot a potential

	ofstream myfile;
	myfile.open("Eigen.txt");

	for(int i=2; i<n_step_width_box; i++){
		//wave_val.push_back(potential_woods(i*h_width +h_width)-centrifug_term(i*h_width +h_width, 0));
		myfile <<i*h_width +h_width<< "   " << 2*U_skyrme_n[i]<< endl;
	//	myfile <<i*h_width +h_width<< "\t" << U_skyrme_p[i] << "\t" << U_skyrme_n[i] << endl;
	}

	myfile.close();



//Plot the proton and neutron central potential
/*
	ofstream myfileN, myfileP;
	myfileN.open("EigenN.txt");
	myfileP.open("EigenP.txt");

	vector<double> Proton_pot;
	vector<double> Neutron_pot;

	for(int i=1; i<n_step_width_box; i++){
		Proton_pot.push_back(potential_woods(i*h_width+h_width));
		Neutron_pot.push_back(potential_woods(i*h_width+h_width)+potential_coulomb(i*h_width+h_width));
		myfileN <<i*h_width +h_width<< "   " << Proton_pot[i] << endl;
		myfileP <<i*h_width +h_width<< "   " << Neutron_pot[i] << endl;
	}

	myfileN.close();
	myfileP.close();
*/

return 0;
}
