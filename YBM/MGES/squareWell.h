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
double numerov(double psi0, double psi1, double E_trial);
double v(double x, double E_trial);
double potential(double x);

// System parameters
extern double a; // well width
extern double hbar2m; // h-bar/2m
extern double h; // step size
extern double epsilon; // precision

// System variables
extern double x; // position in well

#endif
#endif