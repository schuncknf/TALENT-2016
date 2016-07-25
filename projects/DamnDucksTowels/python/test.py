#!/usr/bin/env python

import hffs

basis = hffs.ReducedSpBasis(2.2, 3)
inter = hffs.RawInteraction(basis, 1)
system = hffs.NeutronDrop(basis, inter, 4)
solver = hffs.HartreeFock(system)
solver.initH(0)
for i in range(10):
    solver.run()
    print(solver.cvg)
