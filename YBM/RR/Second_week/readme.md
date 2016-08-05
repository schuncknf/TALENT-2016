
#YBM

3-dimensional spherically symmetric problem for magic nuclei, with Woods-Saxon potential, Coulomb potential and spin-orbit term. 
The program return the eigenvalues and eigenfunctions for the occupied neutron states.
Furthermore, the routine return:
1. The densities of proton and neutron (and the total densities) as function of r.
	(density_neutron.txt, density_proton.txt and density.txt)
2. The potential felt by neutron and proton.
	(EigenN.txt and EigenP.txt)
3. The program also permit to plot all wavefunctions. To have that, you need to uncomment the correpondent part in the main file.


## Usage

First set the values of A, Z and N in the header file. Now it is set for Pb 208.
Just few steps:
1. Make, to compile.
2. Make do, to run the program.
3. Open Gnuplot to obtain the plots. Use: 
	>pl "density_proton.txt" w l, "density_neutron.txt" w l, "density_np.txt" w l
	to see the densities.
	To plot the wavefunctions, if the program has been set to do so, open Gnuplot and write:
	>FILES = system("ls -1 *.dat")
	>plot for [data in FILES] data
	

## Credits

Riccardo Romano
Mattew Shelley
James Cubiss
Acan Dai
