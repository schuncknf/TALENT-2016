#   ____         __                  ___      ___            __          
#  / __/__  ____/ /________ ____    / _ \___ / _/__ ___  ___/ /__ _______
# / _// _ \/ __/ __/ __/ _ `/ _ \  / // / -_) _/ -_) _ \/ _  / -_) __(_-<
#/_/  \___/_/  \__/_/  \_,_/_//_/ /____/\__/_/ \__/_//_/\_,_/\__/_/ /___/

## Coding for Hono(u)r

Welcome to our home in the Repo.

Of particular interest right now is the Hartree-Fock solver in the Hartree-Fock folder.
In the **in** file, you will set the parameters for your HF run.

To see a full set of documentation, visit [here!](http://kylegodbey.com/hf/index.html)

Once you compile and run the code, you'll find the **single particle states**,
**total and kinetic energies**, and **convergence at each iteration** in the file named **out**.

Densities are plotted in the file **densities**, with the order of the columns being *Neutrons, Protons, Total*.

In **plt**, we plot the central parts of the potential for neutrons and protons (in column 1 and 2 respectively)
and then two wavefunctions in the remaining columns. If you want to change which wavefunctions you plot, change it at the bottom of main and recompile.

***Have Fun!***
