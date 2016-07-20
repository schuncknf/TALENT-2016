#include"Esercizio.h"

//Functions included in vector class
/*
	Vettore::Vettore(int a){
		_n=a;
		_v = new double [_n];
	}
	
	Vettore::Vettore(){
		_n=0;
		*_v = 0;
	}
		
	Vettore::~Vettore(){
		delete [] _v;
	}

	void Vettore::Setcomponent(int a, double b){
		_v[a]=b;
	}
	
	void Vettore::Insert(){
		for(int i=0; i<_n;i++)
		cin>>_v[i];
	}
	
	void Vettore::Read(){ 
		for(int i=0; i<_n;i++)
		cout<<_v[i];
	}
	
	double Vettore::Getcomponent(int n){return _v[n];}
*/
//Wave function function

	/*void wave_function(Vettore & x, double energy, double h){
		double appoggio = 2*(1-5*pow(h_bar,2)/12.*v_cost*energy)*x.Getcomponent(1) - (1+pow(h_bar,2)/12.*v_cost*energy)*x.Getcomponent(0)/(1+pow(h_bar,2)/12.*v_cost*energy);		
		x.Setcomponent(0,x.Getcomponent(1));
		x.Setcomponent(1,x.Getcomponent(2));
		x.Setcomponent(2,appoggio);
		//cout<<appoggio<<" ";
	}*/
/*
	Vettore wave_function(Vettore & x, double energy, double h){
		Vettore out(3);
		out.Setcomponent(0,x.Getcomponent(1));
		out.Setcomponent(1,x.Getcomponent(2));
		out.Setcomponent(2,2*(1-5*pow(h_bar,2)/12.*v_cost*energy)*x.Getcomponent(1) - (1+pow(h_bar,2)/12.*v_cost*energy)*x.Getcomponent(0)/(1+pow(h_bar,2)/12.*v_cost*energy));	
		return out;
	}
	*/

		double numerov_algorithm(double energy, double f0, double f_){
			double v = energy/m_factor;
			double a[3];

			a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
			a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
			a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)

			return (a[0] * f0 - a[1] * f_) / a[2];
		}