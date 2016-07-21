#include <stdio.h>
#include <stdlib.h>
#include "gauss-laguerre.h"
// #include "potential.h"
#include "harmon.h"

double fun1 (double x, int n1, int n2, int n3, int n4);

double int1 (double(*pot)(double, int, int, int, int), int n1, int n2, int n3, int n4, double r1, double r2, double m, double w);

double fun2 (double r1, double r2, int n1, int n2, int n3, int n4, double m, double w);

double int2 (double(*int1v)(double(*pot)(double, int, int, int, int),int, int, int, int, double, double, double, double), int n1, int n2, int n3, int n4, double m, double w);

double tbme (double(*pot)(double, double), int n1, int n2, int n3, int n4, double m, double w);
