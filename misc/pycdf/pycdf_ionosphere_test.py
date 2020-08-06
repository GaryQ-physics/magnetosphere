import os
import numpy as np

os.environ["CDF_LIB"] = "/home/gary/theExtractionChamber/CDF_C/cdf38_0-dist/src/lib"
from spacepy import pycdf

#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_130700_000.cdf' # 405.5 kB
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_060700_000.cdf' # 2.7 MB
filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf' # 2.7 MB

cdf = pycdf.CDF(filename)
keys = cdf.keys()

x = cdf['x'][...]
y = cdf['y'][...]
z = cdf['z'][...]
phi = cdf['phi'][...]
lamb = cdf['theta'][...]

x = x[0,:]
y = y[0,:]
z = z[0,:]
phi = phi[0,:]
lamb = lamb[0,:]

x_read = x.reshape((181, 181))
y_read = y.reshape((181, 181))
z_read = z.reshape((181, 181))

deg = np.pi/180.
x_calc = np.outer(np.cos(phi*deg), np.cos(lamb*deg)) # == np.outer(np.cos(lamb*deg), np.cos(phi*deg)).transpose()
y_calc = np.outer(np.sin(phi*deg), np.cos(lamb*deg))
z_calc = np.outer(np.ones(phi.shape), np.sin(lamb*deg))

#print(x_read)
#print(x_calc)
print(np.abs(x_read-x_calc).max())

# x_read==x_calc  ==>  x_read[i,j] == x(lamb[j], phi[i]) where x(lambda, phi) is standard cartesian function of lat, lon (in that order)
#   ==>  var(lamb[j], phi[i]) == cdf['var'][...].reshape((181, 181))[i,j] where variables are considered function var(lambda, phi)
#   ==>  var(lamb[i], phi[j]) == cdf['var'][...].reshape((181, 181)).transpose()[i,j]



def interpolate(lambda_in, phi_in, var): # return var(lambda_in, phi_in) 
    V = cdf[var][...].reshape((181, 181)).transpose()

    #https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
    from scipy.interpolate import RegularGridInterpolator
    ret_interp = RegularGridInterpolator((lamb,phi), V)
    pts = np.column_stack([lambda_in, phi_in])
    return ret_interp(pts)

lambda_in = lamb
phi_in = phi

print(np.column_stack([lambda_in, phi_in]))

x_interp = interpolate(lambda_in, phi_in, 'x')
print(x_interp)
print(x_read)


cdf.close()

if False:
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    fig = plt.figure()
    # https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
    ax = plt.axes(projection='3d') 

    for i in range(181):
        if i%10 != 0: continue
        print(i)
        sol = np.zeros((181, 3))
        for j in range(181):
            sol[j,0] = x_read[i,j]
            sol[j,1] = y_read[i,j]
            sol[j,2] = z_read[i,j]
        #print(sol)
        if i == 0:
            ax.plot3D(sol[:,0], sol[:,1], sol[:,2], 'yellow', lw=4, alpha=1.)
        else:
            ax.plot3D(sol[:,0], sol[:,1], sol[:,2], 'b', lw=2, alpha=0.6)

    for j in range(181):
        if j%10 != 0: continue
        print(j)
        sol = np.zeros((181, 3))
        for i in range(181):
            sol[i,0] = x_read[i,j]
            sol[i,1] = y_read[i,j]
            sol[i,2] = z_read[i,j]
        #print(sol)
        ax.plot3D(sol[:,0], sol[:,1], sol[:,2], 'r', lw=2, alpha=0.6)

    ax.set(xlabel = "X")
    ax.set(ylabel = "Y")
    ax.set(zlabel = "Z")

    L=2
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-L, L)

    plt.show()


