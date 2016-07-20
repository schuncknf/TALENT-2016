#!/usr/bin/env python

import sys
import numpy as np

def main(cfile, hfile, max_degree, step = 10):
    degrees = list(range(0, max_degree + 1, step))
    degrees.remove(0)
    hfile.write("#ifndef QUADRATURES_H\n#define QUADRATURES_H\n\n")
    hfile.write("#define GET_LAG_ROOTS(d, xi, wi) switch(d) { ")
    hfile.write(' '.join(["case %i: xi = lag_p%03i; wi = lag_w%03i; break;" % (d,d,d) for d in degrees]))
    hfile.write(' default: xi = arma::vec(); wi = arma::vec();')
    hfile.write(" }\n")
    
    cfile.write("#include <armadillo>\n\n")
    cfile.write("#include \"quadrature.h\"\n\n")
    cfile.write("// Laguerre quadratures\n")
    for d in degrees:
        p = np.polynomial.Laguerre.basis(d)
        p2 = np.polynomial.Laguerre.basis(d + 1)
        roots = p.roots()
        weights = [ (xi / (d + 1)**2) * (1 / p2(xi))**2 for xi in roots ]
        cfile.write("arma::vec lag_p%03i = { " % d)
        cfile.write(', '.join(["%.15e" % xi for xi in roots]))
        cfile.write(" };\n")
        cfile.write("arma::vec lag_w%03i = { " % d)
        cfile.write(', '.join(["%.15e" % wi for wi in weights]))
        cfile.write(" };\n")
        hfile.write("// Gauss-Laguerre quadrature for d=%i\n" % d)
        hfile.write("extern arma::vec lag_p%03i;\n" % d)
        hfile.write("extern arma::vec lag_w%03i;\n" % d)
    hfile.write("#endif\n")
        
        
        
if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write("usage: %s max [step]\n")
        sys.exit(1)
    hfp = open("quadrature.h", "w")
    cfp = open("quadrature.cpp", "w")
    max_degree = int(sys.argv[1])
    step = int(sys.argv[2]) if len(sys.argv) == 3 else 10
    main(cfp, hfp, max_degree, step)
    hfp.close()
    cfp.close()
