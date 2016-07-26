#ifndef FUNC_H_ 
#define FUNC_H_

#ifndef NS_H_
#define NS_H_

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>

// Functions
double numerov(double x, double psi0, double psi1, double E_trial);
double v(double x, double E_trial);
double potential(double x, int potType);
double checkDLS(double E);

// System parameters
extern double L_box; // box width
extern double L_well; // well width
extern double V0; // well depth
extern int potType; // type of potential well

extern double hbar2m; // h-bar/2m
extern double h; // mesh spacing
extern double epsilon; // precision

#endif
#endif