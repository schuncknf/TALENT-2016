import numpy as npy
from matplotlib import pyplot as plt
from copy import copy, deepcopy
from datetime import datetime as dt
import time

#Import input data
f0 = open(str('VM-scheme_n5_l0.dat')) #reading two body matrix elements.VM-scheme_n5_l2.dat
f1 = open(str('spM_n5_l0.dat'))  #reading orbit number vs quantum number.spM_n5_l2.dat
lines0 = f0.readlines()
lines1 = f1.readlines()
interaction = 1 # 1 - turn on interaction, all other value turns it off
NStates = 12   #108 for n_max = 5, l_max = 2, 12 for n_max = 5, l_max = 0
threshold = 1E-9
#Iteration index
IterMax = 100
partNum = 2
count0 = 0
hw = 10.0 #shorthand for hbar*omega, units in MeV

type = 1    #type 1- neutron, -1- proton, 0- both,
two_iso = type # we only use neutron currently

V0 = {} #dictionary that stores interaction M.E., with key as '(a,b,c,d)'
QNum = {} #dictionary that stores all orbit number in spdata.dat
PNum = {} #Inverse structure as QNum, where keyword becomes value, value becomes keyword
Ind = {}
RInd = {}
#V1 = [[[[0]*N_V]*N_V]*N_V]*N_V

#initialization
E = npy.zeros(NStates)
ePrev = npy.zeros(NStates)
diff = npy.zeros(NStates)
C = npy.zeros([NStates,NStates])
Haml = npy.zeros([NStates,NStates])
H00 = npy.zeros([NStates,NStates])


count = 0
count1 = 0
#Import quantum numbers matrix
for line in lines1:
    count1 += 1
    rs1 = line.split()   #raw string
    try:
        if (count1 < 101):
            OrbNum = int(rs1[2])   #orbit number which corresponds to indices in twobody.dat
            n = int(rs1[3])
            l = int(rs1[4])
            two_j = int(rs1[5])
            two_m = int(rs1[6])
            two_tz = int(rs1[7])
        elif (count1 >= 101):
            str00 = rs1[1].split(":")
            #print str00
            OrbNum = int(str00[1])   #orbit number which corresponds to indices in twobody.dat
            n = int(rs1[2])
            l = int(rs1[3])
            two_j = int(rs1[4])
            two_m = int(rs1[5])
            two_tz = int(rs1[6])
            #Later We'll only use two_tz = 1 for Neutrons
        QNum[(n,l,two_j,two_m,two_tz)] = OrbNum
        PNum[(OrbNum)] = str(n) + ',' + str(l) + ',' + str(two_j) + ',' + str(two_m) + ',' + str(two_tz)
        if (type != 0):
            if (two_tz == two_iso):
                Ind[(OrbNum)] = count
                RInd[(count)] = OrbNum
                count += 1
        else:
            Ind[(OrbNum)] = count
            RInd[(count)] = OrbNum
            count += 1
    except (ValueError, IndexError):
        continue
f1.close()
#print RInd
print "Index Dictionary Complete"


#Import interaction matrix element (M.E.) from tbme.out
for line in lines0:
    rs0 = line.split()   #raw string
    try:
        a = int(rs0[0])
        b = int(rs0[1])
        c = int(rs0[2])
        d = int(rs0[3])
        if ((a in Ind) and (b in Ind) and (c in Ind) and (d in Ind) and (float(rs0[4])!= 0.0)):
            V0[(a,b,c,d)] = float(rs0[4])
    except (ValueError, IndexError):
        continue

f0.close()
#print V1
print "V0 complete"





#define a Kronecker delta function
def KDelta(n1,n3):
    if (n1 == n3):
        return 1.0
    else:
        return 0.0

#density generated eigenvector inputs
def rho(gamma, delta, C):
    rho0 = 0.0
    for N in xrange(0,partNum):
        rho0 += C[gamma][N] * npy.conjugate(C[delta][N])
    return rho0

#alpha=N1,gamma=N2,beta=N3,delta=N4
#units are in MeV, harmonic oscillator basis of hbar*omega = 10 MeV
def spHF(a, b, C):     #need to convert alpha, beta to explicit l,j,n1,n3,
    N1 = RInd[(a)]   #now alpha, beta correspond to Neutron Orbital Number
    N3 = RInd[(b)]    #
    h1 = 0.0
    for i in xrange(0,NStates):
        N2 = RInd[(i)]
        for j in xrange(0,NStates):
            N4 = RInd[(j)]
            if ( ((N1,N2,N3,N4) in V0)):
                h1 += V0[(N1,N2,N3,N4)] * rho(i, j, C)
    
    return h1
 


#use spherical harmonics orthogonality to eliminate some calculation
#either N1 match N3 and N2 match N4; or N1 match N4, N2 match N3.
def spherOrth(N1,N2,N3,N4):
    line1 = PNum[(N1)].split(',')
    l1 = int(line1[1])
    j1 = int(line1[2])
    m1 = int(line1[3])
    line2 = PNum[(N2)].split(',')
    l2 = int(line2[1])
    j2 = int(line2[2])
    m2 = int(line2[3])
    line3 = PNum[(N3)].split(',')
    l3 = int(line3[1])
    j3 = int(line3[2])
    m3 = int(line3[3])
    line4 = PNum[(N4)].split(',')
    l4 = int(line4[1])
    j4 = int(line4[2])
    m4 = int(line4[3])
    if ( (l1 == l3 and j1 == j3 and m1 == m3) and (l2 == l4 and j2 == j4 and m2 == m4) ):
        return 1
    elif ( (l2 == l3 and j2 == j3 and m2 == m3) and (l1 == l4 and j1 == j4 and m1 == m4) ):
        return 1
    else:
        return 0

def eHF(C,partNum):
    eKin = 0
    ePot = 0
    for i in xrange(0,NStates):
        eKin += H00[i][i]*rho(i,i,C)
        for j in xrange(0,NStates):
            for k in xrange(0,NStates):
                for l in xrange(0,NStates):
                    n1 = RInd[(i)]
                    n2 = RInd[(j)]
                    n3 = RInd[(k)]
                    n4 = RInd[(l)]
                    if ( ((n1,n2,n3,n4) in V0)):
                        ePot += 0.5*rho(i,k,C)*rho(j,l,C)*V0[(n1,n2,n3,n4)]

    return eKin+ePot


#making C0 as an Identity, C is eigenvector matrix
for i in xrange(0,NStates):
    C[i][i] = 1.0

#diagonal elements of Hamiltonian
for i in xrange(0,NStates):
    a = RInd[(i)]
    line = PNum[(a)].split(',')
    n = int(line[0])
    l = int(line[1])
    Haml[i][i] = (2*n + l + 1.5) * 10
    H00 = deepcopy(Haml)
print H00

print "Initialization complete, entering HF solver loop..."
#HF loop
Iter = 0
while (Iter < IterMax):
    count3 = 0
    ePrev = deepcopy(E)
    now = dt.now()
    t1 = now.strftime("%M:%S")
    for i in xrange(0,NStates):
        for j in xrange(i,NStates):
            #millis0 = int(round(time.time() * 1000))
            Haml[i][j] += spHF(i, j, C)
            Haml[j][i] = Haml[i][j]
            #millis1 = int(round(time.time() * 1000))
            #delta = millis1 - millis0
            #print '******Time to scan V0 is: ' + str(delta) + ' ms'
    now = dt.now()
    t2 = now.strftime("%M:%S")
    t1_a = t1.split(":")
    t2_a = t2.split(":")
    delta = (int(t2_a[0]) - int (t1_a[0])) * 60 + (int(t2_a[1]) - int (t1_a[1]))

    print '******Time to write Hamiltoninan is: ' + str(delta) + ' seconds'
    print "******Hamiltonian setup complete, begin diagonalization"
    EigValues, EigVectors = npy.linalg.eig(Haml)
    permute = EigValues.argsort()
    EigValues = EigValues[permute]
    EigVectors = EigVectors[:,permute]
    C = deepcopy(EigVectors)
    E = deepcopy(EigValues)
    Haml = deepcopy(H00)
    for i in xrange(0,NStates):
        diff[i] = abs(E[i] - ePrev[i])
    diff.sort()
    Iter += 1
    print 'Iteration: ' + str(Iter) + ', diff = ' + str(diff[NStates-1]) + ' MeV'
    print 'Eigenvalues:'
    print E
    print "Diagonalization complete****************"
    
    if ( diff[NStates-1] < threshold ):
        break

eTot = eHF(C, partNum)
print 'Hartree-Fock Energy = ' + str(eTot) +' MeV\n*\n*\n*'



