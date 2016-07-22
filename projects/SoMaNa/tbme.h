#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"

double fun1 (double, int, int, int, int, double, double, double);

double int1 (double(*pot)(double, int, int, int, int,double,double, double), int, int, int, int, double, double, double, double,double);

double fun2 (double, double, int, int, int, int, double, double,double);

double int2 (double(*int1v)(double,double,int,int,int,int,double,double), int, int, int, int, double, double,double);

double tbme (double(*pot)(double, double), int, int, int, int, double, double);
