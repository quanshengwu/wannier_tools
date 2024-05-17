#!/bin/python3
import numpy as np
import cmath

# The Hamiltonian is 
#     ( M-Bk^2    Delta_0+A*k+  ) 
#     ( Delta_0+A*k_    -M+Bk^2 )
# where k^2=kx^2+ky^2

# phase II, trivial insulator
# A=0, M*B<0

# phase III, trivial insulator with band inversion
# A=0, M*B>0

# phase I, trivial insulator
# Delta_0=0, M*B<0

# phase IV, Chern insulator with band inversion
# Delta_0=0, M*B>0


# from the kp to TB we use sustitution
# k->sin(k)
# k^2->2(1-cos(k))

# Constants
dp = np.float64
pi = np.arctan(1) * 4
zi = 1j

# Lattice constants
M = 2.0
B =-1.0
A = 1.0
Delta_0=  0.0
             

# Number of Wannier functions and R points
num_wann = 4
nrpts = 7

# R coordinates
Irvec = np.zeros((3, nrpts), dtype=int)

# Hamiltonian m,n are band indexes
HmnR = np.zeros((num_wann, num_wann, nrpts), dtype=complex)

# No of degeneracy of R point
ndegen = np.ones(nrpts, dtype=int)

# Initialization of matrices
Irvec[:, :] = 0
HmnR[:, :, :] = 0.0

# 0 0 0
ir = 0
Irvec[:, ir] = [0, 0, 0]
HmnR[0, 0, ir] = M - 4 * B
HmnR[1, 1, ir] = -M + 4 * B
HmnR[2, 2, ir] = M - 4 * B
HmnR[3, 3, ir] = -M + 4 * B
HmnR[0, 1, ir] = Delta_0
HmnR[1, 0, ir] = Delta_0
HmnR[2, 3, ir] = Delta_0
HmnR[3, 2, ir] = Delta_0

# 1 0
ir = 1
Irvec[:, ir] = [1, 0, 0]
HmnR[0, 0, ir] = B
HmnR[1, 1, ir] = -B
HmnR[2, 2, ir] = B
HmnR[3, 3, ir] = -B
HmnR[0, 1, ir] =-0.5*zi*A
HmnR[1, 0, ir] =-0.5*zi*A
HmnR[2, 3, ir] = 0.5*zi*A
HmnR[3, 2, ir] = 0.5*zi*A

# 0 1
ir = 2
Irvec[:, ir] = [0, 1, 0]
HmnR[0, 0, ir] = B
HmnR[1, 1, ir] = -B
HmnR[2, 2, ir] = B
HmnR[3, 3, ir] = -B
HmnR[0, 1, ir]=  -A/2
HmnR[1, 0, ir]=   A/2
HmnR[2, 3, ir]=  -A/2
HmnR[3, 2, ir]=   A/2


# -1 0
ir = 3
Irvec[:, ir] = [-1, 0, 0]
HmnR[0, 0, ir] = B
HmnR[1, 1, ir] = -B
HmnR[2, 2, ir] = B
HmnR[3, 3, ir] = -B
HmnR[0, 1, ir]= 0.5*zi*A
HmnR[1, 0, ir]= 0.5*zi*A
HmnR[2, 3, ir]=-0.5*zi*A
HmnR[3, 2, ir]=-0.5*zi*A


# 0 -1
ir = 4
Irvec[:, ir] = [0, -1, 0]
HmnR[0, 0, ir] = B
HmnR[1, 1, ir] = -B
HmnR[2, 2, ir] = B
HmnR[3, 3, ir] = -B
HmnR[0, 1, ir]=   A/2
HmnR[1, 0, ir]=  -A/2
HmnR[2, 3, ir]=   A/2
HmnR[3, 2, ir]=  -A/2

nrpts= ir+1
# Writing to a file
with open('HalfBHZ_hr.dat', 'w') as file:
    file.write('2-band half of BHZ model\n')
    file.write('2\n')
    file.write(f'{nrpts}\n')
    file.write(' '.join(f'{x:5d}' for x in ndegen) + '\n')
    for ir in range(nrpts):
        for i in range(2):
            for j in range(2):
                file.write(f"{Irvec[0, ir]:5d}{Irvec[1, ir]:5d}{Irvec[2, ir]:5d}{i+1:5d}{j+1:5d} {HmnR[i, j, ir].real:16.8f} {HmnR[i, j, ir].imag:16.8f}\n")

