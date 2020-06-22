"""
Demonstrates how kameleon can work with odeint

odeint is a function in scipy.integrate which requires a callable function and then solves the differential equation associated with it

a callable function that accesses kameleon can be used if we pass the kameleon and interpolator objects as constant arguments

to show this we solve the differential equation for a pendulum with extra force term given by the pressure 'p' in kameleon
"""

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

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

def Compute(plot_force_function = False):
    # open kameleon ---------------
    kameleon = ccmc.Kameleon()
    kameleon.open(filename)
    print(filename, "Opened " + filename)
    interpolator = kameleon.createNewInterpolator()
    #-----------------------------

    A=np.linspace(-30.,30.,300)
    S=np.zeros((300,))
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
    if plot_force_function:
        plt.plot(A, S, 'r', label='F(theta)')
        plt.xlabel('theta')
    else:
        plt.plot(t, sol[:, 0], 'b', label='theta(t)')
        plt.plot(t, sol[:, 1], 'g', label='omega(t)')
        plt.xlabel('t')

    plt.legend(loc='best')
    plt.grid()
    plt.show()

    # close kameleon ---------------------
    kameleon.close()
    print("Closed " + filename)
    #-------------------------------

    return [sol, S]

OUT=Compute()
