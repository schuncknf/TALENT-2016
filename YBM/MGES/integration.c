#include "squareWell.h"


double a = 6.0;
double hbar2m = 20.75;
double h = 0.1;
double epsilon = 0.001;


double numerov(double psi0, double psi1, double x, double E_trial){
    double a1 = 2.0*(1.0-(pow(h,2)*5.0/12.0)*v(x, E_trial));
    double a2 = (1.0+(pow(h,2)/12.0))*v(x-h, E_trial);
    double a3 = (1.0+(pow(h,2)/12.0))*v(x+h, E_trial);
    
    return (a1*psi1 - a2*psi0) / a3;
}

double v(double x, double E_trial){
    return (E_trial-potential(x)) / hbar2m;
}

double potential(double x){
    return 0.0;
}