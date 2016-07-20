#include "squareWell.h"
using namespace std;

double numerov(double psi0, double psi1, double E_trial){
    double a1 = 2.0*(1.0 - h*h*5.0/12.0*v(x,E_trial));
    double a2 = (1.0 + h*h/12.0*v(x-h,E_trial));
    double a3 = (1.0 + h*h/12.0*v(x+h,E_trial));
    
    return (a1*psi1 - a2*psi0) / a3;
}


double v(double x, double E_trial){
    return (E_trial-potential(x)) / hbar2m;
}


double potential(double x){
    return 0.0;
}