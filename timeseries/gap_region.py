import numpy as np
from numba import njit
from swmf_file_reader.read_swmf_files import read_all, interpolate
import datetime
import pandas as pd

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from magnetometers import GetMagnetometerLocation
from units_and_constants import phys
import cxtransform as cx

print(cx.transform([1.,0.,0.], (2019,9,2,6,30,0), 'GSM', 'MAG'))#!!!!!!!! DOESNBT WORK UNLESS UNCOMMENTED
print(cx.transform([1.,0.,0.], (2019,9,2,6,30,0), 'GSM', 'SM'))
print(cx.transform([0.,1.,0.], (2019,9,2,6,30,0), 'GSM', 'SM'))
print(cx.transform([0.,0.,1.], (2019,9,2,6,30,0), 'GSM', 'SM'))

rCurrents = 1.8
rIonosphere = 1.01725 # rEarth + iono_height
GM_csys = 'GSM'

@njit
def sph_to_xyz(r, theta, phi):
    Xyz_D = np.empty(3)
    Xyz_D[0] = r*np.cos(phi)*np.sin(theta)
    Xyz_D[1] = r*np.sin(phi)*np.sin(theta)
    Xyz_D[2] = r*np.cos(theta)
    return Xyz_D

@njit
def map_along_dipole_lines(Xyz_D, rMap):
    # Xyz_D and returned XyzMap_D in SMG (SM) coordinates
    # Solution of the vector potential equation (proportional to (x^2+y^2)/r^3)
    # so sqrt(xMap^2+yMap^2)/sqrt(x^2+y^2) = sqrt(rMap^3/r^3)
    iHemisphere = int(np.sign(Xyz_D[2]))
    XyzMap_D = np.empty(3, dtype=np.float32)

    r = np.linalg.norm(Xyz_D)
    XyRatio = np.sqrt(rMap/r)**3 # ratio of input and mapped X and Y components
    XyzMap_D[0:2] = XyRatio*Xyz_D[0:2]
    XyMap2 = XyzMap_D[0]**2 + XyzMap_D[1]**2

    if rMap**2 < XyMap2:
       # The point does not map to the given radius
       iHemisphere = 0
       # Put mapped point to the magnetic equator
       XyzMap_D[0:2] = (rMap/np.sqrt(Xyz_D[0]**2 + Xyz_D[1]**2))*Xyz_D[0:2]
       XyzMap_D[2] = 0
    else:
       XyzMap_D[2] = iHemisphere*np.sqrt(rMap**2 - XyMap2)

    return XyzMap_D, iHemisphere

@njit
def get_dipole_field(xyz):
    # Xyz_D and returned b_D in SMG (SM) coordinates
    b = np.empty(3)
    r = np.linalg.norm(xyz)
    DipoleStrength = 3.12e+4 #"dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field
    Term1      = DipoleStrength*xyz[2]*3/r**2
    b[0:2] = Term1*xyz[0:2]/r**3
    b[2]    = (Term1*xyz[2]-DipoleStrength)/r**3
    return b

@njit
def _jit_initialize(nTheta, nPhi, nR):
    xyz_Currents_V = np.empty((nTheta, nPhi, 3))
    dSurface_V = np.empty((nTheta, nPhi))

    dTheta = np.pi    / (nTheta-1)
    dPhi   = 2.*np.pi / nPhi
    dR     = (rCurrents - rIonosphere) / nR
    for iTheta in range(nTheta):
        Theta = iTheta * dTheta
        # the area of the triangle formed by pole and the latitude segment at Theta=dTheta/2
        # is approximately dTheta/4*dTheta/2, so sin(theta) replaced with dTheta/8.
        SinTheta = max(np.sin(Theta), dTheta/8.)
        dSurface = rCurrents**2*SinTheta*dTheta*dPhi
        for iPhi in range(nPhi):
            Phi = iPhi * dPhi

            dSurface_V[iTheta, iPhi] = dSurface
            xyz_Currents_V[iTheta, iPhi] = sph_to_xyz(rCurrents, Theta, Phi)

    return xyz_Currents_V, dSurface_V

@njit
def _jit_loop(xyz_Currents_V, dSurface_V, b_Currents_V, j_Currents_V, x0, nTheta, nPhi, nR):
    dB_fac                 = np.zeros(3)
    dB_mhd_SurfaceIntegral = np.zeros(3)

    dTheta = np.pi    / (nTheta-1)
    dPhi   = 2.*np.pi / nPhi
    dR     = (rCurrents - rIonosphere) / nR
    for iTheta in range(nTheta):
        for iPhi in range(nPhi):
            dSurface = dSurface_V[iTheta,iPhi]
            xyz_Currents = xyz_Currents_V[iTheta,iPhi,:]
            b0_Currents = get_dipole_field(xyz_Currents)
            b_Currents = b_Currents_V[iTheta,iPhi,:]
            j_Currents = j_Currents_V[iTheta,iPhi,:]

            Unit_xyz_Currents = xyz_Currents / rCurrents
            Unit_b_Currents = b_Currents / np.linalg.norm(b_Currents)
            _Fac_rCurrents = np.sum(Unit_b_Currents*j_Currents)#!!!!!!!!!!!
            Fac_term = _Fac_rCurrents * np.sum(Unit_b_Currents*Unit_xyz_Currents)

            #####################
            Br   = np.sum(Unit_xyz_Currents*b_Currents)
            Bt = np.cross(Unit_xyz_Currents, b_Currents)
            InvDist2_D = dSurface*(xyz_Currents - x0)/(4*np.pi*np.sqrt(np.sum((xyz_Currents - x0)**2))**3)
            dB_mhd_SurfaceIntegral[:] = dB_mhd_SurfaceIntegral[:] + Br*InvDist2_D + np.cross(Bt, InvDist2_D)
            #####################

            for k in range(nR):
                R = rCurrents - dR*(k+0.5)

                xyz_Map, iHemisphere = map_along_dipole_lines(xyz_Currents, R)
                b0_Map = get_dipole_field(xyz_Map)

                Unit_xyz_Map = xyz_Map / R
                Unit_b0_Map = b0_Map / np.linalg.norm(b0_Map)

                # The volume element is proportional to 1/Br. The sign
                # should be preserved (not yet!!!),
                # because the sign is also there in the radial
                # component of the field aligned current: Br/B*FAC.
                # In the end j_D = b_D/Br*[(Br/B)*(j.B)]_rcurr  !!!!!!!!!!!! NOT DIMENSIONALLY CONSISTENT

                dVol_FACcoords = dSurface * (dR/np.sum(Unit_xyz_Map*Unit_b0_Map))
                J_fac = Fac_term * Unit_b0_Map

                dB_fac[:] = dB_fac[:] + dVol_FACcoords* \
                  np.cross(J_fac, x0-xyz_Map)/(4*np.pi*(np.linalg.norm(xyz_Map-x0))**3)

    return dB_fac, dB_mhd_SurfaceIntegral


def _vectorized_gap_region_integrals(cache, time, x0, nTheta,nPhi,nR, gap_csys='SM'):
    # x0 and returned dB_fac, dB_mhd_SurfaceIntegral are in cartesian gap_csys coordinates (default SMG)
    b_Currents_V = np.empty((nTheta,nPhi,3), dtype=np.float32)
    j_Currents_V = np.empty((nTheta,nPhi,3), dtype=np.float32)

    xyz_Currents_V, dSurface_V = _jit_initialize(nTheta, nPhi, nR)
    trans_xyz_Currents_V = cx.transform(xyz_Currents_V.reshape((nTheta*nPhi,3)), time, gap_csys, GM_csys)
    b_Currents_V[:,:,0] = interpolate(None, trans_xyz_Currents_V, var='bx', cache=cache).reshape((nTheta,nPhi))
    b_Currents_V[:,:,1] = interpolate(None, trans_xyz_Currents_V, var='by', cache=cache).reshape((nTheta,nPhi))
    b_Currents_V[:,:,2] = interpolate(None, trans_xyz_Currents_V, var='bz', cache=cache).reshape((nTheta,nPhi))
    j_Currents_V[:,:,0] = interpolate(None, trans_xyz_Currents_V, var='jx', cache=cache).reshape((nTheta,nPhi))
    j_Currents_V[:,:,1] = interpolate(None, trans_xyz_Currents_V, var='jy', cache=cache).reshape((nTheta,nPhi))
    j_Currents_V[:,:,2] = interpolate(None, trans_xyz_Currents_V, var='jz', cache=cache).reshape((nTheta,nPhi))
    b_Currents_V[:,:,:] = cx.transform(b_Currents_V.reshape((nTheta*nPhi,3)), time,GM_csys,gap_csys).reshape((nTheta,nPhi,3))
    j_Currents_V[:,:,:] = cx.transform(j_Currents_V.reshape((nTheta*nPhi,3)), time,GM_csys,gap_csys).reshape((nTheta,nPhi,3))

    dB_fac, dB_mhd_SurfaceIntegral = _jit_loop(xyz_Currents_V, dSurface_V, b_Currents_V, j_Currents_V, x0, nTheta, nPhi, nR)
    dB_fac = (phys['mu0']*phys['muA']/phys['m']**2) * dB_fac
    return dB_fac, dB_mhd_SurfaceIntegral

def slice_gap_region(run, time, obs_point, nTheta=181,nPhi=180,nR=30, gap_csys='SM', cache=None):
    assert(gap_csys=='SM')
    if cache is None:
        cache = read_all(util.time2CDFfilename(run,time)[:-8])

    if isinstance(obs_point,str):
        obs_point_str = obs_point
        if obs_point == "origin":
            x0 = np.zeros(3)
        else:
            x0 = GetMagnetometerLocation(obs_point_str, time, 'SM', 'car')
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])
        x0 = cx.transform(obs_point, time, 'GSM', 'SM')
    x0 = np.array(x0)

    outname_fac = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'B_fac_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s.npy'%(obs_point_str)

    outname_surf = conf[run+'_derived'] + 'timeseries/slices/' \
        + 'B_mhd_SurfaceIntegral_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
        + '_obs_point=%s.npy'%(obs_point_str)


    dB_fac, dB_mhd_SurfaceIntegral = _vectorized_gap_region_integrals(cache, time, x0, nTheta,nPhi,nR, gap_csys=gap_csys)
#    print(dB_fac)
#    print(dB_mhd_SurfaceIntegral)
    np.save(outname_fac, dB_fac)
    np.save(outname_surf, dB_mhd_SurfaceIntegral)

def stitch_gap_region(run, times, obs_point, nTheta=181,nPhi=180,nR=30, gap_csys='SM'):
    assert(gap_csys=='SM')

    if isinstance(obs_point,str):
        obs_point_str = obs_point
    else:
        obs_point_str = '[%f,%f,%f]'%(obs_point[0],obs_point[1],obs_point[2])

    df_name = conf[run+'_derived']+'timeseries/df_gap_region' \
            + '_obs_point=%s.pkl'%(obs_point_str)

    dtimes = []
    slice_arrays = []
    for time in times:
        dtimes.append(datetime.datetime(*time[:6]))
        outname_fac = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'B_fac_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_obs_point=%s.npy'%(obs_point_str)
        outname_surf = conf[run+'_derived'] + 'timeseries/slices/' \
            + 'B_mhd_SurfaceIntegral_%.2d%.2d%.2dT%.2d%.2d%.2d'%util.tpad(time, length=6) \
            + '_obs_point=%s.npy'%(obs_point_str)

        dB_fac = np.load(outname_fac)
        dB_mhd_SurfaceIntegral = np.load(outname_surf)

        slice_arrays.append(np.concatenate([dB_fac,dB_mhd_SurfaceIntegral]))

    columns = ('B_fac_x','B_fac_y','B_fac_z', 'B_mhd_SurfaceIntegral_x', 'B_mhd_SurfaceIntegral_y', 'B_mhd_SurfaceIntegral_z')
    df = pd.DataFrame(data=slice_arrays, columns = columns,
                        index=dtimes)
    df.to_pickle(df_name)

if __name__ == '__main__':
    run = 'DIPTSUR2'
    #time = (2019,9,2,6,30,0)
    #obs_point = np.array([0.,0.,0.])
    #slice_gap_region(run,time,obs_point)
    #stitch_gap_region(run,list(util.get_available_slices(run)[1]), "colaba")
