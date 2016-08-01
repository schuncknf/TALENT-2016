#!/usr/bin/env python

import hffs
import config
import numpy as np


def read_TBME(inter, f_TBME):
    basis = []
    nMax = 3
    lMax = 1
    jMin = 0
    jMax = 5
    for n in range(1, nMax + 1):
        for l in range(0, lMax + 1):
            for j2 in range(int(2 * jMin + 1), int(2 * jMax), 2):
                for tz in range(-1, 1):
                    for mj2 in range(-j2, j2 + 1, 2):
                        basis.append((n,l,j2 + 1,mj2,tz))
                        print (n,l,j2,mj2,tz)
    for l in fp:
        if i in (0, 1):
            continue
        
        

if __name__ == '__main__':
    if config.lMax == 0:
        basisType = "ReducedSpBasis"
        basis = hffs.ReducedSpBasis(config.omega, config.nMax)
    else:
        basisType = "SpBasis"
        basis = hffs.SpBasis(config.omega, config.nMax, config.lMax, config.mMax)

    pi = config.interaction.split(':')
    if len(pi) == 1:
        inter_input = "direct"
        inter_type, = pi
    elif len(pi) == 2:
        inter_input, inter_type = pi
    else:
        raise ValueError("Bad syntax")

    if inter_input == "direct":
        if inter_type == "MinnesotaS0":
            if basisType == "ReducedSpBasis":
                inter = hffs.MinnesotaS0(basis, 60)
            else:
                raise ValueError("Bad basis type \"%s\" for \"%s\" interaction" % (basisType, config.interaction))
        else:
            raise ValueError("Unknown interaction \"%s\"" % (config.interaction))
    elif inter_input == "file":
        pass
        #fp = open(inter_type, "r")
        # TODO
    else:
        raise ValueError("Unknown input \"%s\"" % (input_type))
        
    if config.system == "NeutronDrop":
        system = hffs.NeutronDrop(basis, inter, config.nb_neutron, config.omega, 150)
    elif config.system == "Nucleus":
        system = hffs.Nucleus(basis, inter, config.nb_neutron, config.nb_proton)
    else:
        raise ValueError("Unknown system \"%s\"" % (config.system))

    if config.solver in ("Hartree-Fock", "HF"):
        solver = hffs.HartreeFock(system)
    elif config.solver in ("Hartree-Fock-Bogoliubov", "Hartree-Fock-Bogo", "HFB"):
        solver = hffs.HartreeFockBogoliubov(system)
    else:
        raise ValueError("Unknown solver \"%s\"" % (config.solver))

    solver.initH(0)
    cvg = 1000
    i = 0
    np.set_printoptions(linewidth = 1000, suppress = True)
    print("%03i %s" % (i,solver.info()))
    while cvg > config.convergence:
        solver.run()
        system.calcH()
        cvg = solver.cvg
        i += 1
        rho = system.getR(0,0)
        kinetic = system.getKinetic(0)
        gamma = system.getGamma(0)
        delta = system.getDelta(0)
        #print(solver.indivEnergies)
        ene = np.trace(kinetic.dot(system.getR(0,0)) + 0.5 * gamma.dot(system.getR(0,0)))
        print("%03i %s, e=%.7f" % (i,solver.info(), ene))
