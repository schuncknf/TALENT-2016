#!/usr/bin/env python

import hffs

basis = hffs.ReducedSpBasis(2.2, 3)
inter = hffs.RawInteraction(basis, 1)
system = hffs.NeutronDrop(basis, inter, 4)
solver = hffs.HartreeFock(system)
for i in range(10):
    solver.calc()
    print(solver.cvg)
