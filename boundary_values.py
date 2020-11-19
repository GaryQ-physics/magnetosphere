import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import util
from probe import probe
import cxtransform as cx
from units_and_constants import phys

#sys.path.append('/home/gary/magnetovis/pkg/magnetovis/')
#from vtk_export import vtk_export


run = 'DIPTSUR2'
pntlist = 'native_random_sampled'

if run == 'DIPTSUR2':
    time = (2019,9,2,6,30,0,0)
    rCurrents = 1.8
if run == 'IMP10_RUN_SAMPLE':
    time = (2019,9,2,7,0,0,0)
    rCurrents = 1.7
if run == 'TESTANALYTIC':
    time = (2000,1,1,0,10,0,0)
    rCurrents = 0.

direct = conf[run+'_derived'] + 'regions/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
direct = direct + pntlist + '/'
points = np.loadtxt(direct + pntlist + '_points.txt')
if debug: print(points)

direct = conf[run+'_derived'] + 'boundaries/%.2d%.2d%.2dT%.2d%.2d%.2d/'%util.tpad(time, length=6)
if not os.path.exists(direct):
    os.makedirs(direct)

import time as tm
t0 = tm.time()

deg = np.pi/180.

R = rCurrents
Nlat = 2*180 + 1
Nlon = 2*360 + 1

lat_ax = np.linspace(-90., 90., Nlat)
lon_ax = np.linspace(0., 360., Nlon)

dtheta = deg*180./(Nlat-1)
dphi = deg*360./(Nlon-1)

lon, lat = np.meshgrid(lon_ax, lat_ax)
P_sph = np.column_stack([R*np.ones(Nlat*Nlon), lat.flatten(), lon.flatten()])
P = cx.transform(P_sph, time, 'MAG', 'GSM', ctype_in='sph', ctype_out='car')

J = probe(util.time2CDFfilename(run, time), P, var=['jx','jy','jz'], library='kameleon')
#J *= phys['mu0']*(phys['muA']/(phys['m']**2))

print(np.all(np.sqrt(np.einsum('ij,ij->i',P,P)) == R))
J_perp = np.einsum('ij,ij->i', J, P)/R
J_norm = np.sqrt(np.einsum('ij,ij->i', J, J))

#outfname = conf[run+'_derived'] + 'J_perp_boundary_R_%s.vtk'%(R)

arr = np.column_stack([P, J_perp, J_norm])
np.savetxt('J_perp.txt', arr)


dA = np.zeros((Nlat,Nlon))
dA[0,0] = np.pi*(dtheta/2.)**2
dA[Nlat-1,0] = np.pi*(dtheta/2.)**2
toinsert = np.cos(lat_ax[1:Nlat-1]*deg)*dtheta*dphi
toinsert = np.tile(toinsert,(Nlon-1,1)).transpose()
dA[1:Nlat-1, 0:Nlon-1] = toinsert

dA *= R**2

#J_perp=J_perp.reshape((Nlat,Nlon))
dA = dA.flatten()

if False:
    I = np.sum(J_norm*dA)
    print(I)

    I = np.sum(np.abs(J_perp)*dA)
    print(I)

    I = np.sum(J_perp*dA)
    print(I)

I = np.nan*np.empty(points.shape)
for i in range(points.shape[0]):
    x = points[i,:]
    rvect = P - x
    rdist = np.einsum('ij,ij->i',rvect, rvect)
    Fx = J_perp*rvect[:,0]/(rdist**3)
    Fy = J_perp*rvect[:,1]/(rdist**3)
    Fz = J_perp*rvect[:,2]/(rdist**3)
    I[i,0] = np.sum(Fx*dA)
    I[i,1] = np.sum(Fy*dA)
    I[i,2] = np.sum(Fz*dA)

np.savetxt(direct + 'bound.txt', np.column_stack([points, I]))


'''
vtk_export(outfname, P,
                    dataset = 'STRUCTURED_GRID',
                    connectivity = (Nlat,Nlon,1),
                    point_data = J_perp,
                    texture = 'SCALARS',
                    point_data_name = 'J_perp',
                    )
'''
