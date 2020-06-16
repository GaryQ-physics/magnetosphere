# vtk_functions

import sys
import os
import numpy as np
import meshio

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from units_and_constants import phys
import biot_savart as bs

def writevtk(fname, p_np, cust_cells):
    c_np = [('hexahedron', cust_cells.astype(int))]
    fvals_dict= {'sample_scalars': fvals}
    cust_mesh = meshio.Mesh(p_np, c_np, point_data=fvals_dict)
    print("Writing " + fname + 'with meshio')
    meshio.vtk.write(fname, cust_mesh, binary=False)
    print("Wrote " + fname)

if True:
    fname = conf["run_path_derived"] + 'vtk_functions_test.vtk'

    dx = 0.05*phys['R_e']
    dy = 0.05
    dz = 0.1
    Rad = 0.5
    xlims = [-1.1*Rad, 1.1*Rad]
    ylims = [-1.1*Rad, 1.1*Rad]
    zlims = [-1.1*Rad, 1.1*Rad]

    ret = bs.make_grid(xlims, ylims, zlims, dx, dy, dz)
    p_np = ret[0]
    fvals = 0.1*p_np[:,0]*p_np[:,1]*p_np[:,2]
    cell_inds = ret[4]
    ind_p_np = ret[5]
    Nx = ret[1]
    Ny = ret[2]
    Nz = ret[3]

    shifts = np.array( [[0, 0, 0],
                        [1, 0, 0],
                        [1, 1, 0],
                        [0, 1, 0],
                        [0, 0, 1],
                        [1, 0, 1],
                        [1, 1, 1],
                        [0, 1, 1]] )

    cust_cells = np.nan*np.empty(((Nx-1)*(Ny-1)*(Nz-1), 8))
    test_cells = np.nan*np.empty(((Nx-1)*(Ny-1)*(Nz-1), 8))
    for k in range((Nx-1)*(Ny-1)*(Nz-1)):
        for n in range(8):
            tup = tuple(cell_inds[k,:] + shifts[n,:])
            #vtup = tuple(cell_inds[k,:] + shifts[n,:])
            cust_cells[k,n] = np.where(np.all([ind_p_np[:,0]==tup[0],ind_p_np[:,1]==tup[1],ind_p_np[:,2]==tup[2]],axis=0))[0][0]
            #test_cells[k,n] = np.where(np.all([p_np[:,0]==vtup[0],p_np[:,1]==vtup[1],p_np[:,2]==vtup[2]],axis=0))[0][0]

    print(cust_cells)
    print(cust_cells.shape)

    writevtk(fname, p_np, cust_cells)
