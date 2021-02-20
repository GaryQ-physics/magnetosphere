import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../../')
import jitFORTRAN

with open('subroutine.f', 'r') as include:
    script = ''.join(include.readlines())

#print(script)

meadV_F = jitFORTRAN.Fortran_Subroutine(script, 'MEADV')

x = np.arange(10, dtype=np.float64)+1.
y = np.arange(10, dtype=np.float64)+1.
z = np.arange(10, dtype=np.float64)+1.

bx = np.nan*np.empty((10,), dtype=np.float64)
by = np.nan*np.empty((10,), dtype=np.float64)
bz = np.nan*np.empty((10,), dtype=np.float64)

meadV_F.execute(x, y, z, 1, bx, by, bz, 10)

print(bx, by, bz)
