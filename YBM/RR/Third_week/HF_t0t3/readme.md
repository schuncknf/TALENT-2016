
#YBM

HF calculation with Skyrme t0t3 functional. The code is written to run with Oxygen 16. <br />
The routine works with an initial Woods-Saxon potential. 
The program return the eigenvalues and eigenfunctions for the occupied neutron states.

## Usage

Just few steps:
1. Make, to compile.
2. Make do, to run the program.
To add the Coulomb potential, you have to uncomment the corresponded part just after the loop starts (line 197) and uncomment the terms when the Numerov algorithm is called (line 308, 331, 341).
The total energy doesn't include the coulomb energy.
	

## Credits

Riccardo Romano
Mattew Shelley
James Cubiss
Acan Dai
