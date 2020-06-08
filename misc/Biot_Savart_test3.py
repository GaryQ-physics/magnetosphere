# Biot_Savart_test3

import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

# units
hr = 1.
muA = 1.
R_e = 1.
nT = 1.

deg = (np.pi/180.)
amin = deg/60.
minn = hr/60.
s = minn/60.

A = 1e+6*muA
Tesla = 1e+9*nT
m = R_e/6.3781e+6  # R_e == 6.3781e+6*m  (note m < 10^-6 Re)
#kg = Tesla*A*s**2

mu0_SI = 1.25663706212e-6 # almost 4*pi*1e-7
if False:
    mu0 = mu0_SI*kg*m/((s**2)*(A**2))
else:
    mu0 = 1.970237314e-10*nT*R_e/muA


Npole = np.array([0., 0., 1.])
Rad = 0.5 # considered in units R_e
j_mag = 1. # considered in units muA/m^2
def J_kunits(x, y, z):
    X = np.array([x, y, z])
    if X[0]**2 + X[1]**2 < Rad**2:
        j = j_mag*np.array([0,0,1]) #in units of muA/m^2
    else:
        j = np.zeros((3,))
    return j

Itot = np.pi*(Rad*R_e)**2*(j_mag*muA/m**2) # pi * (0.5*6.3781e+6 m)**2 * (1. muA/m^2) = 3.19501226e+13 muA = 3.19501226e+10 Amps

print('Itot*1e-13= ', Itot*1e-13, ' ?=? 3.19501226')

def interpolate(variable, x, y, z):
    j = J_kunits(x, y, z)
    if variable == 'jx':
        return j[0]
    if variable == 'jy':
        return j[1]
    if variable == 'jz':
        return j[2]

def ex_data_full(kam,interp, variable, x, y, z, X0, Npole, V_char = 1.):
    #if np.sqrt(x**2+y**2+z**2)<1e-4: return 0.
    # Get data from file, interpolate to point
    if('dB' in variable):
        #if np.sqrt(x**2+y**2+z**2)<1.: return 0.
        if np.dot(X0,X0)<1e-8: return 0.
        J = np.array([ex_data_full(kam, interp, 'jx', x, y, z, X0, Npole), 
                      ex_data_full(kam, interp, 'jy', x, y, z, X0, Npole), 
                      ex_data_full(kam, interp, 'jz', x, y, z, X0, Npole)])
        J = J*(muA/m**2)
        R = X0 - np.array([x, y, z])
        #R = R*R_e
        #dB_dV = (mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
        #dBnT = dB_dV*V_char/(nT)
        dBnT = V_char*(mu0/(4*np.pi))*np.cross(J, R)/(np.linalg.norm(R)**3)
        if(variable=='dB'):
            return np.linalg.norm(dBnT)
        a2 = np.cross(Npole, X0)
        if(variable=='dB_EW'):
            return np.dot(dBnT, a2)/np.linalg.norm(a2) # east west direction (east positive)
        a1 = np.cross(X0, a2)
        if(variable=='dB_NS'):
            return np.dot(dBnT, a1)/np.linalg.norm(a1) # north south direction (north positive)
    data = interpolate(variable, x, y, z)
    return data


def B_EW(X0, Npole, mult=1, length=10.):
    dx = 0.05/mult
    dy = 0.05/mult
    dz = 0.1/mult
    X = np.arange(-Rad - 0.1, Rad + 0.1 + dx, dx)
    Nx = X.size # 24
    Y = np.arange(-Rad - 0.1, Rad + 0.1 + dy, dy)
    Ny = Y.size
    Z = np.arange(-0.5*length, 0.5*length + dz, dz)
    Nz = Z.size # 100
    B2, B3, B1 = np.meshgrid(Y, Z, X)
    #B1, B2, B3 = np.meshgrid(X, Y, Z)

    B1 = B1.flatten(order='C')
    B2 = B2.flatten(order='C')
    B3 = B3.flatten(order='C')
    Bgrid = np.column_stack((B1, B2, B3))


    total = 0.

    #X0 = np.array(X0)

    #Aa = (np.nan)*np.empty((B1.size, ))
    for l in range(B1.size):
        #Aa[l] = ex_data_full(False, False, 'dB_EW', B[l, 0], B[l, 1], B[l, 2], X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3
        total = total + ex_data_full(False, False, 'dB_EW', Bgrid[l, 0], Bgrid[l, 1], Bgrid[l, 2], X0, Npole, V_char = dx*dy*dz) # dx*dy*dz*R_e**3
        #print('total=',total)
    return total


def B_analytic(X0):
    X0 = np.array(X0)
    r = np.linalg.norm(X0)

    r_SI = r*6.3781e+6
    Rad_SI = Rad*6.3781e+6
    J_SI = j_mag*1e-6
    I_SI = J_SI*np.pi*Rad_SI**2 # 3.19501226e+7

    Itot = np.pi*(Rad*R_e)**2*(j_mag*muA/m**2) # pi * (0.5*6.3781e+6 m)**2 * (1. muA/m^2) = 3.19501226e+13 muA = 3.19501226e+7 Amps
    #print('Itot= ', Itot)
    #print('I_SI= ', I_SI)
    #print(' ?=? 3.19501226')

    B_SI = mu0_SI*I_SI/(2*np.pi*r_SI) # 1.33582614*1e-6
    return B_SI*1e+9  # norm of analytic B in nanotesla


if True:
    print B_EW([0., 0.75, 0.], Npole) # note for this input EW direction is parallel to expected analytic B
    print B_analytic([0., 0.75, 0.])

'''
x0 = np.linspace(-1., 1., 2)
y0 = np.linspace(-1., 1., 2)

X0g, Y0g = np.meshgrid(x0, y0)

Bn = (np.nan)*np.empty(X0g.shape)
Bn_a = (np.nan)*np.empty(X0g.shape)

for i in range(X0g.shape[0]):
    for j in range(X0g.shape[1]):
        Bn[i,j] = B_EW([X0g[i,j], Y0g[i,j], 0.])**2 + B_NS([X0g[i,j], Y0g[i,j], 0.])**2
        #Bn_a[i,j] = np.dot( B_analytic([X0g[i,j], Y0g[i,j], 0.]), B_analytic([X0g[i,j], Y0g[i,j], 0.]) )

print Bn
print Bn_a
print X0g
print Y0g
'''

