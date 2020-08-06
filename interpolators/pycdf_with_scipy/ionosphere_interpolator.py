import os
import numpy as np

os.environ["CDF_LIB"] = "/home/gary/theExtractionChamber/CDF_C/cdf38_0-dist/src/lib"
from spacepy import pycdf
from scipy.interpolate import RegularGridInterpolator


def interpolate(lamb, phi, var, filename): # return var(lambda_in, phi_in) 
    cdf = pycdf.CDF(filename)
    #keys = cdf.keys()

    lamb_ax = cdf['theta'][...][0,:]
    phi_ax = cdf['phi'][...][0,:]

    # x_read==x_calc  ==>  x_read[i,j] == x(lamb[j], phi[i]) where x(lambda, phi) is standard cartesian function of lat, lon (in that order)
    #   ==>  var(lamb[j], phi[i]) == cdf['var'][...].reshape((181, 181))[i,j] where variables are considered function var(lambda, phi)
    #   ==>  var(lamb[i], phi[j]) == cdf['var'][...].reshape((181, 181)).transpose()[i,j]

    V = cdf[var][...].reshape((181, 181)).transpose()

    cdf.close()

    #https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
    ret_interp = RegularGridInterpolator((lamb_ax,phi_ax), V)
    pts = np.column_stack([lamb, phi])

    return ret_interp(pts)

