
#YBM

3-dimensional spherically symmetric problem for magic nuclei, with Woods-Saxon potential, Coulomb potential and spin-orbit term. <br />
The program return the eigenvalues and eigenfunctions for the occupied neutron states. <br />
Furthermore, the routine return: <br />
1. The densities of proton and neutron (and the total densities) as function of r. <br />
	(density_neutron.txt, density_proton.txt and density.txt) <br />
2. The potential felt by neutron and proton. <br />
	(EigenN.txt and EigenP.txt) <br />
3. The program also permit to plot all wavefunctions. To have that, you need to uncomment the correpondent part in the main file. 


## Usage

First set the values of A, Z and N in the header file. Now it is set for Pb 208. <br />
Just few steps: <br />
1. Make, to compile. <br />
2. Make do, to run the program. <br />
3. Open Gnuplot to obtain the plots. Use: <br />
	>pl "density_proton.txt" w l, "density_neutron.txt" w l, "density_np.txt" w l <br />
	to see the densities. <br />
	To plot the wavefunctions, if the program has been set to do so, open Gnuplot and write: <br />
	>FILES = system("ls -1 *.dat") <br />
	>plot for [data in FILES] data
	

## Credits

Riccardo Romano <br />
Mattew Shelley <br />
James Cubiss <br />
Acan Dai  <br />
