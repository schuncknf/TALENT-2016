nox = True

# Basis configuration
omega = 0.0506779468
nMax = 5
lMax = 2
mMax = 0

interaction = "file:./minnesota-VM-scheme.dat"
#interaction = "file:/run/shm/minnesota-VM-scheme.dat"
system = "NeutronDrop"
nb_neutron = 8
#nb_proton = 5

solver = "Hartree-Fock"
convergence = 1.0e-13
