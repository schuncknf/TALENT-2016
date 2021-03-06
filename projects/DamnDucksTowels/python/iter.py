#!/usr/bin/env python

import shf
import config
import numpy as np
import time

if __name__ == '__main__':
    if not config.nox:
        logofp = open("team-logo.txt", "r")
        ss = logofp.read()
        print(ss)
        logofp.close()
        logofp = open("small-code-logo.txt", "r")
        ss = ""
        print("+" * 85 + "\n")
        for l in logofp:
            vl = len(l)
            ss += l
        print(ss)
        logofp.close()
        print("+" * 85 + "\n")
    ############################################
    t0 = time.time()
    ############################################
    if config.lMax == 0:
        basisType = "ReducedSpBasis"
        basis = shf.ReducedSpBasis(config.omega, config.nMax)
    else:
        basisType = "FullSpBasis"
        basis = shf.FullSpBasis(config.omega, config.nMax, config.lMax)

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
                inter = shf.MinnesotaS0(basis, 60)
            else:
                raise ValueError("Bad basis type \"%s\" for \"%s\" interaction" % (basisType, config.interaction))
        else:
            raise ValueError("Unknown interaction \"%s\"" % (config.interaction))
    elif inter_input == "file":
        inter = shf.MinnesotaRaw(basis, config.nb_neutron, inter_type)
        #fp = open(inter_type, "r")
        # TODO
    else:
        raise ValueError("Unknown input \"%s\"" % (input_type))
        
    if config.system == "NeutronDrop":
        system = shf.NeutronDrop(basis, inter, config.nb_neutron, config.omega, 150)
    elif config.system == "Nucleus":
        system = shf.Nucleus(basis, inter, config.nb_neutron, config.nb_proton)
    else:
        raise ValueError("Unknown system \"%s\"" % (config.system))

    if config.solver in ("Hartree-Fock", "HF"):
        solver = shf.HartreeFock(system)
    elif config.solver in ("Hartree-Fock-Bogoliubov", "Hartree-Fock-Bogo", "HFB"):
        solver = shf.HartreeFockBogoliubov(system)
    else:
        raise ValueError("Unknown solver \"%s\"" % (config.solver))

    solver.initH(0)
    cvg = 1000
    i = 0
    np.set_printoptions(linewidth = 1000, suppress = True)
    ############################################
    t1 = time.time()
    ############################################
    if config.nox:
        while cvg > config.convergence:
            solver.run()
            system.calcH()
            cvg = solver.cvg
            i += 1
            #print(kinetic + gamma.dot(rho))
        rho = system.getR(0,0)
        kinetic = system.getKinetic(0)
        gamma = system.getGamma(0)
        #delta = system.getDelta(0)
        #print(solver.indivEnergies)
        ene = np.trace(kinetic.dot(system.getR(0,0)) + 0.5 * gamma.dot(system.getR(0,0)))
        print("%03i %s, e=%.7f" % (i,solver.info(), ene))
    else:
        print("%03i %s" % (i,solver.info()))
        while cvg > config.convergence:
            solver.run()
            system.calcH()
            cvg = solver.cvg
            i += 1
            rho = system.getR(0,0)
            kinetic = system.getKinetic(0)
            gamma = system.getGamma(0)
            #delta = system.getDelta(0)
            #print(solver.indivEnergies)
            ene = np.trace(kinetic.dot(system.getR(0,0)) + 0.5 * gamma.dot(system.getR(0,0)))
            print("%03i %s, e=%.7f" % (i,solver.info(), ene))
            #print(kinetic + gamma.dot(rho))
        
        
    ############################################
    t2 = time.time()
    ############################################
    
    print(t1 - t0)
    print(t2 - t1)
