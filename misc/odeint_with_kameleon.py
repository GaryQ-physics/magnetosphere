# kameleon_test

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()

import _CCMC as ccmc
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Inputs for file and magnetic field model
year = 2003
day = 20
month = 11
hours = 7
minutes = 0
seconds = 0

filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % (year,month,day,hours,minutes,seconds) + '.out.cdf'



def ex_data(kam,interp,variable, x,y,z):
    # Get data from file, interpolate to point
    #kam.loadVariable(variable)
    data = interp.interpolate(variable, x, y, z)
    if( x**2 + y**2 + z**2 >=1.):
        return data
    else:
        return 0.

def pend(y, t, b, c, kam, interp,):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta) - ex_data(kam,interp,'p',omega,omega,omega)]
    return dydt

def Compute():
    # open kameleon ---------------
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()
    #-----------------------------

    A=np.linspace(1.,2.,100)
    S=np.zeros((100,))
    #print A

    kameleon.loadVariable('p')

    for k in range(A.size):
        #print 'k=',k
        S[k]=ex_data(kameleon,interpolator,'p',A[k],A[k],A[k])

    b = 0.25
    c = 5.0
    y0 = [np.pi - 0.1, 0.0]
    t = np.linspace(0, 10, 101)
    sol = odeint(pend, y0, t, args=(b, c, kameleon, interpolator))
    plt.plot(t, sol[:, 0], 'b', label='theta(t)')
    plt.plot(t, sol[:, 1], 'g', label='omega(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------

    return S


def writevtk():
    out_fname=conf["m_path"] + 'magnetosphere/data/' + 'test_out'+'.vtk'
    sol=Compute()
    f = open(out_fname,'w')
    print('writing ' + out_fname)
    f.write('# vtk DataFile Version 3.0\n')

    for k in range(sol.size):
        f.write('%e\n'%(sol[k]))

    f.close()
    print('closed ' + out_fname)

#writevtk()
OUT=Compute()
