import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../../')
import jitFORTRAN

with open('outer_sub.f', 'r') as include:
    script = ''.join(include.readlines())

#print(script)

meadV_F = jitFORTRAN.Fortran_Subroutine(script, 'MEADV', include='sub')

N=10

x = np.arange(N, dtype=np.float64)+1.
y = np.arange(N, dtype=np.float64)+1.
z = np.arange(N, dtype=np.float64)+1.

bx = np.nan*np.empty((N,), dtype=np.float64)
by = np.nan*np.empty((N,), dtype=np.float64)
bz = np.nan*np.empty((N,), dtype=np.float64)

meadV_F.execute(x, y, z, 1, bx, by, bz, N)

print(bx, by, bz)
