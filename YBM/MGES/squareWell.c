#include "squareWell.h"


int main(){
    
    double E_trial = 1.0;
    double x=0.0;
    
    double psi0 = 0.0;
    double psi1 = 1.0;
    
    std::cout << x << psi0 << std::endl;
    
    // Loop from left to right of well 
    for(int i=1;i<=a/h;++i){
        double psi_temp = psi1;
        double psi1 = numerov(psi0,psi1,x,E_trial);
        psi0 = psi_temp;
        x += h;
        std::cout << x << "\t" << psi0 << std::endl;
    }
    return 0;
}

