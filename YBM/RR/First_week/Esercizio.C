#include "Esercizio.h"

int main(){

	double a = 35;
	double h = 0.1;
	double n_step = a/h;

	Vettore x(3);
	x.Setcomponent(0,0.);
	x.Setcomponent(1,1.);

	//for(int j=1; j<10; j++){
		double energy = 0.1671;
		//Vettore x_new(3);
		for(int i=0; i<n_step; i++){
			x = wave_function(x, energy, h);
			//x = x_new;
			cout<<x.Getcomponent(0)<<endl;
		}
		//cout<<x.Getcomponent(0)<<endl;
	//}
return 0;
}
