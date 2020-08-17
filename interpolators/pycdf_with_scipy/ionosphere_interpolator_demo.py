import numpy as np
import ionosphere_interpolator as ii


filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'
lamb = np.array([20, 10, 5, 1, 0., -1, -5, -10, -20])
phi = np.array([0., 0, 0, 0, 0, 0, 0, 0, 0]) + 10.

ret = ii.interpolate(lamb, phi, 'jx', filename)

print(ret)
#print(ret==None)

