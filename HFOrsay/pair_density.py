# -*- coding: utf-8 -*-
#!/usr/bin/python
#!/usr/bin/env python

import os
import sys
import subprocess
#cmd = './main < /dev/null'
cmd = './run'
FNULL = open(os.devnull, 'w')
if not os.path.exists('pairing_density'):
    os.makedirs('pairing_density')
pair_i =  0.5
while pair_i < 3.0:
    pair_i = pair_i + 0.02  
    input_file = open('hforsay.dat','w')
    input_file.write("n0base   =   5              ! Number of osc. shells\n")
    input_file.write("ngauss   =  120             ! Points for the laguerre-Mesh\n")
    input_file.write("npart   =    8              ! Number of particules\n")
    input_file.write("maxit    =   30             ! Maximum iterations number\n")
    input_file.write("homega   =    10.00\n")
    input_file.write("bcs      =    1              ! Plug in pairing\n")
    input_file.write("gpair    =     %.4f\n" % pair_i)
    input_file.write("plot     =    1              ! Ploting Routines\n")
    input_file.write("iflag    =    1              ! From external tbme.bin\n")
    input_file.write("i3d      =    1              ! From external tbme.bin\n")
    input_file.close()
    subprocess.call(['./run'])       
    if os.path.exists('density.dat'):
        os.rename('density.dat','pairing_density/density_%.6f' % pair_i +".dat")
        os.remove('hforsay.dat')
#            if os.path.exists('s010_010.hel'):
#                os.remove('s010_010.hel')
