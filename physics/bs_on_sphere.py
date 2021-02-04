import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
deg = np.pi/180.

from probe import probe
import biot_savart as bs
import cxtransform as cx

time = [2006,12,5,0,0,0]

x0 = cx.MAGtoGSM(np.array([1., 45.03, 45.]), time, 'sph', 'car')

Nlat = 2*180 + 1
Nlon = 2*360 + 1

lat_ax = np.linspace(-90., 90., Nlat)
lon_ax = np.linspace(0., 360., Nlon)

dtheta = deg*180./(Nlat-1)
dphi = deg*360./(Nlon-1)

#lat, lon = np.meshgrid(lat_ax, lon_ax)
lon, lat = np.meshgrid(lon_ax, lat_ax)

P = np.column_stack([lat.flatten(), lon.flatten()])

filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'

units = 1. #!!!!!!!!!!!!!!
J = probe(filename, P, var=['jx','jy','jz'], library='pycdf')*units  # assuming x,y,z comp of J are in GSM
X = cx.MAGtoGSM(np.column_stack([np.ones(Nlat*Nlon), P]), time, 'sph', 'car')

print(lat.shape)
print(lat.flatten().shape)
print(P.shape)
print(X.shape)

u_vect = np.array([1., 0., 0.])
dB = bs.deltaB('dB', x0, X, J, V_char = 1.)

f = np.einsum('ij,j', dB, u_vect).reshape((Nlat,Nlon))
#f = probe(filename, P, var='jy', library='pycdf')
#f = np.ones((Nlat, Nlon))

I = 0
for i in range(Nlat):
    if i == 0 or i == Nlat-1:
        dA = np.pi*(dtheta/2.)**2
        I += dA*f[i, 0]

    else:
        for j in range(Nlon - 1):
            dA = np.cos(lat_ax[i]*deg)*dtheta*dphi
            I += dA*f[i,j]


dA = np.zeros((Nlat,Nlon))
dA[0,0] = np.pi*(dtheta/2.)**2
dA[Nlat-1,0] = np.pi*(dtheta/2.)**2
toinsert = np.cos(lat_ax[1:Nlat-1]*deg)*dtheta*dphi
toinsert = np.tile(toinsert,(Nlon-1,1)).transpose()
dA[1:Nlat-1, 0:Nlon-1] = toinsert

I2 = np.sum(dA*f)

I3vect = bs.deltaB('deltaB', x0, X, J, V_char = dA.flatten())

I3 = np.dot(I3vect, u_vect)

'''
dA_ = np.zeros((Nlat,Nlon))
for i in range(Nlat):
    if i == 0 or i == Nlat-1:
        dA_[i,0] = np.pi*(dtheta/2.)**2
    else:
        for j in range(Nlon - 1):
            dA_[i,j] = np.cos(lat_ax[i]*deg)*dtheta*dphi

print(np.all(dA_ == dA)) --> True
'''



print(I) # 9.090509756157932e-06
print(I2)
print(I3)
print(4*np.pi)


filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'
time = [2006,12,5,0,0,0]
x0 = cx.MAGtoGSM(np.array([1., 45.03, 45.]), time, 'sph', 'car')

def ionosphere_biot_savart(filename, time, x0, Nlat=2*180 + 1, Nlon=2*360 + 1):
    lat_ax = np.linspace(-90., 90., Nlat)
    lon_ax = np.linspace(0., 360., Nlon)
    dtheta = deg*180./(Nlat-1)
    dphi = deg*360./(Nlon-1)

    lon, lat = np.meshgrid(lon_ax, lat_ax)
    P = np.column_stack([lat.flatten(), lon.flatten()])

    units = 1. #!!!!!!!!!!!!!!
    J = probe(filename, P, var=['jx','jy','jz'], library='pycdf')*units  # assuming x,y,z comp of J are in GSM
    X = cx.MAGtoGSM(np.column_stack([np.ones(Nlat*Nlon), P]), time, 'sph', 'car')

    dA = np.zeros((Nlat,Nlon))
    dA[0,0] = np.pi*(dtheta/2.)**2
    dA[Nlat-1,0] = np.pi*(dtheta/2.)**2
    toinsert = np.cos(lat_ax[1:Nlat-1]*deg)*dtheta*dphi
    toinsert = np.tile(toinsert,(Nlon-1,1)).transpose()
    dA[1:Nlat-1, 0:Nlon-1] = toinsert

    ret = bs.deltaB('deltaB', x0, X, J, V_char = dA.flatten())

