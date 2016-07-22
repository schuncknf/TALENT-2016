#include "Nucleus.h"

Nucleus::Nucleus (Basis & basis, int nbNeut, int nbProt, arma::field<arma::mat> & TBME) 
    : System("Nucleus", basis, arma::ivec {nbNeut,nbProt}, std::vector<std::string> {"Neutrons","Protons"}, TBME)
{}

void Nucleus::calcH0(arma::field<arma::mat> & H, int type) {
    (void) H;
    (void) type;
}

void Nucleus::calcH(arma::field<arma::mat> & H, arma::field<arma::mat> & RG) {
    (void) H;
    (void) RG;
}

void Nucleus::calcKineticField(arma::field<arma::mat> & RG){
    (void) RG;
}

Nucleus::~Nucleus(){}
