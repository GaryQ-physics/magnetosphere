import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import odeint

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf
#import _CCMC as ccmc

import pos_sun as ps
import b_field_lines_write as bfl


########################################################################
def ex_data(kam, interp, variable, x, y, z):
    """Load data from file, interpolate to point"""
    if (x**2 + y**2 + z**2 < 1.):
        return 0
    data = interp.interpolate(variable, x, y, z)
    return data

def dXds(X, s, kam, interp, var, sign):
    """Derivative function for field line ODE

    dx/ds = Fx(x,y,z)/Fm
    dy/ds = Fy(x,y,z)/Fm
    dz/ds = Fz(x,y,z)/Fm
    
    X = [x, y, z]
    F = [Fx, Fy, Fz]
    Fm = sqrt(Fx**2 + Fy**2 + Fz**2)
    s = arclength

    F is magnetic field for          var = 'b'
    F is current density field for   var = 'j'
    """
    
    F = np.array([ex_data(kam, interp, var + 'x', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'y', X[0], X[1], X[2]), 
                  ex_data(kam, interp, var + 'z', X[0], X[1], X[2])])
    Fm = np.sqrt(np.dot(F, F))
    if 1e-9 < Fm < 1e+7:
        return (sign/Fm)*F
    else:
        return [0., 0., 0.]

def Compute(Event, Nb, debug=False):
    #Event = [year, month, day, hours, minutes, seconds, miliseconds, MLONdeg, MLATdeg]
    time = Event[0:7]
    MLON = Event[7]
    MLAT = Event[8]
    T = tuple(time)

    filename = conf["run_path"] + '3d__var_3_e' \
                + '%04d%02d%02d-%02d%02d%02d-%03d' % T + '.out.cdf'

    R = 1.01
    eps = 3.
    IC = []
    for i in range(Nb+1):
        if i==0:
            delta = 0.
        elif i>Nb/2:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        else:
            #delta=(i-Nb/2)*eps
            delta = (-i)*eps
        IC.append(ps.MAGtoGSM([R, MLAT-delta, MLON], time[0:6], 'sph', 'car'))

    if debug:
        print(IC)

    # Trace field lines
    s_grid = np.linspace(0., 200., 2000)

    if debug:
        print(s_grid)
        print(s_grid.size)

    solns = (np.nan)*np.empty((s_grid.size, 3, len(IC)))
    for i in range(len(IC)):
        kameleon = ccmc.Kameleon()
        if debug:
            print("Opening " + filename)
        kameleon.open(filename)
        if debug:
            print("Opened " + filename)
        interpolator = kameleon.createNewInterpolator()
        kameleon.loadVariable('bx')
        kameleon.loadVariable('by')
        kameleon.loadVariable('bz')
        sol = odeint(dXds, IC[i], s_grid, args = (kameleon, interpolator, 'b', -1))
        kameleon.close()
        solns[:, :, i] = sol

    if debug:
        print('IC', IC)

    # restrict the field lines to stop when reaching 1*R_E from the origin
    solns_restr = [] # initialize list of np_arrays, one for each restricted field line
    for i in range(Nb+1):  # loop over field lines
        # define condition on the field line points
        #tr = np.logical_and(solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 >=1.,solns[:,0,i]**2+solns[:,1,i]**2+solns[:,2,i]**2 <= 224.**2+128.**2+128.**2)
        # create the arrays of the restricted field line componentwise
        went_out = 0
        end_val = solns.shape[0]-1
        if debug:
            print(end_val)
        for k in range(solns.shape[0]):
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 > 1.1**2: went_out = 1
            if solns[k, 0, i]**2 + solns[k, 1, i]**2 + solns[k, 2, i]**2 < 1.1**2 and went_out==1:
                end_val = k
                break
        if debug:
            print(end_val)
        tr = np.arange(end_val+1)
        solx = solns[:, 0, i]
        if debug:
            print(solx)
            print(tr)
        solx = solx[tr]

        if debug:
            print(solx)
        soly = solns[:, 1, i]
        soly = soly[tr]
        solz = solns[:, 2, i]
        solz = solz[tr]

        # reasemble and add to the list
        sol = np.column_stack([solx, soly, solz])
        solns_restr.append(sol)
        if debug and i==Nb+1: print(solns[:, :, i])
        if debug and i==Nb+1: print(sol)

    if debug:
        print(solns_restr)
    return solns_restr
########################################################################


print(sys.version)
print(sys.path)
print(os.environ['PYTHONPATH'].split(os.pathsep))

fig = plt.figure()
ax = plt.axes(projection='3d') #https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

solns_original = Compute([2003, 11, 20 , 7, 0, 0, 0, 176.00, 57.50], 6)
solns_withGrid = bfl.Compute([2003, 11, 20, 7, 0, 57.50, 176.00], 6)

print(len(solns_withGrid) == len(solns_original))

sol0 = solns_original[0]
sol1 = solns_withGrid[0]

for i in range(len(solns_withGrid)):
    sol0 = solns_original[i]
    sol1 = solns_withGrid[i]

    min_len = min(sol0.shape[0], sol1.shape[0])
    print(np.abs(sol1[:min_len,:] - sol0[:min_len,:]).max())

    print(sol0.shape)
    print(sol1.shape)
    print(sol0[0,:])
    print(sol1[0,:])

    ax.plot3D(sol0[:,0], sol0[:,1], sol0[:,2], 'red', lw=2)
    ax.plot3D(sol1[:,0], sol1[:,1], sol1[:,2], 'blue', lw=4, alpha=0.6)

ax.set(xlabel = "$X/R_E$ (GSM)")
ax.set(ylabel = "$Y/R_E$ (GSM)")
ax.set(zlabel = "$Z/R_E$ (GSM)")

L=20
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ax.set_zlim(-L, L)

plt.show()
