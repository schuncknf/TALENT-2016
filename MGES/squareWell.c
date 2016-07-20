#include "squareWell.h"

using namespace std;


double a = 6.0;
double hbar2m = 20.75;
double h = 0.1;
double epsilon = 0.001;

double E_trial = 10.;
double x = 0.0;
double psi0 = 0.0;
double psi1 = 1.0;


int main(){
    
    // Output
    std::ofstream ofile;
    ofile.open("output.dat");

    
    ofile << x << "\t" << psi0 << endl;
    
    // Loop from left to right of well 
    for(int i=1;i<=a/h;++i){
        
        // Store last calculated wavefn value in 'temp' variable
        double psi_temp = psi1;
        // Calculate next wavefn value
        psi1 = numerov(psi0,psi1,E_trial);
        // Set previous value 'psi0'
        psi0 = psi_temp;
        // Step across well
        x += h;
        
        ofile << x << "\t" << psi0 << endl;
        
    }
    
    ofile.close();
    
    return 0;
}

