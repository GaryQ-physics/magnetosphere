import os
import numpy as np

os.environ["CDF_LIB"] = "/home/gary/theExtractionChamber/CDF_C/cdf38_0-dist/src/lib"
from spacepy import pycdf
from scipy.interpolate import RegularGridInterpolator


def interpolate(lamb, phi, var, filename): # return var(lambda_in, phi_in) 
    cdf = pycdf.CDF(filename, readonly=True)
    #keys = cdf.keys()

    lamb_ax = cdf['theta'][...][0,:]
    phi_ax = cdf['phi'][...][0,:]

    lamb_ax[0] = -90.
    lamb_ax[-1] = 90.  # before its like 89.999...something but needs to
                       # be 90 otherwise when calling interpolator about
                       # this 89.99some, doesnt work. Needed fore ex in bs_on_sphere.py

    # x_read==x_calc  ==>  x_read[i,j] == x(lamb[j], phi[i]) where x(lambda, phi) is standard cartesian function of lat, lon (in that order)
    #   ==>  var(lamb[j], phi[i]) == cdf['var'][...].reshape((181, 181))[i,j] where variables are considered function var(lambda, phi)
    #   ==>  var(lamb[i], phi[j]) == cdf['var'][...].reshape((181, 181)).transpose()[i,j]

    V = cdf[var][...].reshape((181, 181)).transpose()

    cdf.close()

    pts = np.column_stack([lamb, phi])
    if pts = np.column_stack(np.meshgrid(lamb_ax, phi_ax)): #!!!!!!!!!!
        print("on grid points")
        return pts

    #https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
    ret_interp = RegularGridInterpolator((lamb_ax,phi_ax), V)
    return ret_interp(pts)

