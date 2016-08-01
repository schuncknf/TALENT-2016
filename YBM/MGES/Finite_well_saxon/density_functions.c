double Coulomb(double r){
    double Rp = 1.0; // Proton radius
    double e = 1.6e-19; // Charge constant
    
    // Z is defined in header file

    if(r <= Rp){
        return Z*e*e/(2.*Rp)*(3. - pow(r/Rp,2));
    }
    else{
        return Z*e*e/r
    }

}

double<vector> singleOrbitDensity(double j, double<vector> efn){
    // Returns radial density function for single eigenfunction 'efn'
    
    double<vector> density(efn.size());
        
    for(int i=0;density.size();++i){
        density[i] = (2*j + 1)/(4*pi*pow(i*h,2) * pow(efn[i],2);
    }
    
    return density;
}

double totalMatterDensity(double<vector> density){
    // Integrate total radial density 'density' to calculate A
    double A = 0.0;
    
    for(i=1;density.size();++i){
        A += (density[i]-density[i-1])/h;
    }
    
    return A;
}