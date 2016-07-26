#include"Esercizio.h"

	double numerov_algorithm(double energy, double f0, double f_, double x, int potType, int l, double s){
		double v = (energy-potential(x, potType, l, s))/m_factor;
		double a[3];
		a[0] = 2. * (1. - 5./12. * v * pow(h_width,2));	// Coeff. for f(x)
		a[1] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x-h)
		a[2] = 1. * (1. + 1./12. * v * pow(h_width,2));	// Coeff. for f(x+h)
		return (a[0] * f0 - a[1] * f_) / a[2];
	}


	double potential(double x, int potType, int l, double s){
        double     WS = coeffV * (1 / ( 1 + expFactor));
        int         A = nProton + nNeutron;
        double      R = r0 * pow(A,1./3);
        double coeffV = -51 + 33*((nNeutron-nProton)/A);
        double expFactor = exp( (x-(R/2))/a_comp );
        
		if(potType == 0){
            return 0;
        }

        if(potType == 1){
            double half_box = (width-width)/2.;
            if(x<half_box)
                return V0;
            else if(x>=half_box&&x<=6+half_box)
                return 0;
            else if(x>6+half_box)
                return V0;
        }
        
        else if(potType == 2){
            return WS;
        }
        
        else if(potType == 3){            
            double dfdr = -1.*pow(1. + expFactor,-2) * expFactor;
            
            return WS*Vls*l*s*pow(r0,2)*dfdr/x;
        }
        
	}