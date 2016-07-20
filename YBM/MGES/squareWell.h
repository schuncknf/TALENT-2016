#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>


// Functions
double numerov(double psi0, double psi1, double x, double E_trial);
double v(double x, double E_trial);
double potential(double x);

// Globals
extern double a; // well width
extern double hbar2m; // h-bar/2m
extern double h; // step size
extern double epsilon; // precision