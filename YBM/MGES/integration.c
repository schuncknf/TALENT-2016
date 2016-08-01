#include "squareWell.h"
using namespace std;


double numerov(double x, double psi0, double psi1, double E_trial){
    // Implements Numerov algorithm
    double a1 = 2.0*(1.0 - h*h*5.0/12.0*v(x,E_trial));
    double a2 = (1.0 + h*h/12.0*v(x-h,E_trial));
    double a3 = (1.0 + h*h/12.0*v(x+h,E_trial));
    
    return (a1*psi1 - a2*psi0) / a3;
}


double v(double x, double E_trial){
    return (E_trial-potential(x,potType)) / hbar2m;
}


double potential(double x, int potType){
    if(potType == 0){ // Infinite square well - use positive energies
        if(x > (L_box-L_well)/2.0 && x < (L_box+L_well)/2.0){
            return 0.0;
        }
        else{
            return 1e10; // Effectively infinite for an energy in [0-100 MeV]
        }
    }
    else if(potType == 1){ // Finite square well
        if(x >= (L_box-L_well)/2.0 && x <= (L_box+L_well)/2.0){
            return -V0;
        }
        else{
            return 0.0;
        }
    }
    else{
        cout << "Choose well type between 0 and 1" << endl;
        return 0;
    }
}


double checkDLS(double E_trial){
    
    // Output
    ofstream DLSfile;
    DLSfile.open("DLS.dat");

    double x=0.0;
    double psi0 = 0.0;
    double psi1 = 1.0e-5;
    
    // Integrate from LHS of box to RHS of well
    for(int i=1;i<=(L_box+L_well)/2.0/h;++i){
        
        DLSfile << x << "\t" << psi0 << endl;
        
        // Store last calculated wavefn value in 'temp' variable
        double psi_temp = psi1;
        // Calculate next wavefn value
        psi1 = numerov(x,psi0,psi1,E_trial);
        // Set previous value 'psi0'
        psi0 = psi_temp;
        // Step across well
        x += h;
    }
    
    // Store final LHS value and derivative
    double LHSval = psi1;
    double LHSderiv = (psi1-psi0)/h;
    
    // Integrate from RHS of box to RHS of well
    x=L_box;
    psi0 = 0.0;
    psi1 = 1.0e-5;
    
    h = -h;
    
    for(int i=L_box/-h;i>=(L_box+L_well)/2.0/-h;--i){
        
        DLSfile << x << "\t" << psi0 << endl;    
        
        // Store last calculated wavefn value in 'temp' variable
        double psi_temp = psi1;
        // Calculate next wavefn value
        psi1 = numerov(x,psi0,psi1,E_trial);
        // Set previous value 'psi0'
        psi0 = psi_temp;
        // Step across well
        x += h;
    }    
    
    // Store final RHS value and derivative
    double RHSval = psi1;
    double RHSderiv = (psi0-psi1)/-h;
    
    DLSfile.close();
    
    cout << RHSderiv/LHSderiv - RHSval/LHSval << endl;
    
    // Return DLS condition
    return RHSderiv/LHSderiv - RHSval/LHSval;
    
}