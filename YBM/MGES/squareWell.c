#include "squareWell.h"

using namespace std;


double hbar2m = 20.75;
double h = 0.01;
double epsilon = 0.001;

double L_well = 35.0;
double E_trial = -2.5;
double L_box = 59.0;
double V0 = 20.0;
// 0 > inf; 1 > finite square well; 
int potType = 1;

double x = 0.0;
double psi0 = 0.0;
double psi1 = 1.0e-5;


int main(){
    
    // Output
    ofstream ofile;
    ofile.open("output.dat");

    ofile << x << "\t" << psi0 << "\t" << potential(x, potType) << endl;
    
    // Loop from left to right of well
    for(int i=1;i<=L_box/h;++i){
        
        // Store last calculated wavefn value in 'temp' variable
        double psi_temp = psi1;
        // Calculate next wavefn value
        psi1 = numerov(x,psi0,psi1,E_trial);
        // Set previous value 'psi0'
        psi0 = psi_temp;
        // Step across well
        x += h;
        
        ofile << x << "\t" << psi0 << "\t" << potential(x, potType) << endl;
        
    }
    
    ofile << x << "\t" << psi1 << endl;
    
    ofile.close();
    
    cout << checkDLS(E_trial) << endl;
    
    return 0;
}

