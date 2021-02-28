# Depreciated

import os
import sys
sys.path.append('/home/gary/magnetosphere/physics/')
from make_grid import make_axes, make_grid
import numpy as np
import spacepy.pybats.bats as bats
import read_swmf_files as rswmf

#https://stackoverflow.com/questions/2706605/sorting-a-2d-numpy-array-by-multiple-axes
def isReordering(arr1, arr2, debug=False):
    arr1 = np.array(arr1); arr2 = np.array(arr2)
    if arr1.shape != arr2.shape: return False

    sort1 = np.lexsort((arr1[:,0],arr1[:,1],arr1[:,2]))  
    sort2 = np.lexsort((arr2[:,0],arr2[:,1],arr2[:,2]))
    undo1 = np.argsort(sort1)
    undo2 = np.argsort(sort2)

    if debug:
        print('#######')
        print(np.all(np.arange(5896192) == sort1[undo1]))
        print(np.all(np.arange(5896192) == sort2[undo2]))
        print(np.all(arr1 == (arr1[sort1,:])[undo1,:]))
        print(np.all(arr2 == (arr2[sort2,:])[undo2,:]))
        print(np.all(arr2 == (arr1[sort1,:])[undo2,:]))
        print(np.all( arr2 == (arr1[sort1[undo2],:]) ))
        print('#######')

    return np.all(arr1[sort1,:] == arr2[sort2,:])

def Reorder(arr1, arr2):
    ''' Reorder(arr1, arr2), if isReordering(arr1, arr2) is true,
            returns array of indices 'ind' s.t  arr1[ind,:] == arr2
    arr1 = np.array(arr1); arr2 = np.array(arr2)
    '''
    arr1 = np.array(arr1); arr2 = np.array(arr2)
    if arr1.shape != arr2.shape: return False

    sort1 = np.lexsort((arr1[:,0],arr1[:,1],arr1[:,2]))  
    sort2 = np.lexsort((arr2[:,0],arr2[:,1],arr2[:,2]))
    #undo1 = np.argsort(sort1)
    undo2 = np.argsort(sort2)

    if not np.all(arr2 == (arr1[sort1[undo2],:])):
        raise RuntimeError ('arrays inputed arent reordering of each other')

    return sort1[undo2]


def generate_indices(fname, testOUTvsVTU=False, testOUTvsRECON=False, save=True):
    '''
    testOUTvsVTU=True only works with python3
    '''

    data = bats.Bats2d(fname + ".out")
    points = np.column_stack([data['x'],data['y'],data['z']])

    if testOUTvsVTU:
        import pyvista as pv
        mesh = pv.read(fname + ".vtu")
        assert(np.all(8*np.arange(mesh.offset.size) == mesh.offset))
        assert(np.all(mesh.celltypes == 12))

        print('points from .out file:')
        print(points)
        print('points from .vtu file:')
        print(mesh.points)
        print('is points from .out and .vtu the same points in the same order:')
        print(np.all(points == mesh.points))
        print('is points from .out and .vtu the same points in the potentially different order:')
        print(isReordering(points, mesh.points))

    npts = points.shape[0]
    allgrids = []
    gridspacings = []
    counter = 0
    dTREE = rswmf.read_tree_file(fname)
    blockused_ = rswmf.get_blocks_used(dTREE)
    nI, nJ, nK = dTREE['nI'], dTREE['nJ'], dTREE['nK']
    for iNode in blockused_:
        xlims, ylims, zlims, gridspacing = rswmf.get_physical_dimensions(iNode, dTREE, returnCenters=True)

        assert(xlims[1]-xlims[0] == (nI-1)*gridspacing)
        assert(ylims[1]-ylims[0] == (nJ-1)*gridspacing)
        assert(zlims[1]-zlims[0] == (nK-1)*gridspacing)

        #### way 1 ####
        #grid = np.empty((nI*nJ*nK,3))
        #cnt=0
        #for i in range(nI):
        #    for j in range(nJ):
        #        for k in range(nK):
        #            grid[cnt, 0] = xlims[0]+gridspacing*i
        #            grid[cnt, 1] = ylims[0]+gridspacing*j
        #            grid[cnt, 2] = zlims[0]+gridspacing*k
        #            cnt += 1

        #### way 2 ####
        #grid = np.empty((nI,nJ,nK,3))
        #for i in range(nI):
        #    for j in range(nJ):
        #        for k in range(nK):
        #            grid[i,j,k, 0] = xlims[0]+gridspacing*i
        #            grid[i,j,k, 1] = ylims[0]+gridspacing*j
        #            grid[i,j,k, 2] = zlims[0]+gridspacing*k
        #grid = grid.reshape((nI*nJ*nK,3))

        #### way 3 ####
        grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                        ylims[0]:ylims[1]+gridspacing:gridspacing,
                        zlims[0]:zlims[1]+gridspacing:gridspacing ]
        grid = np.array(grid.reshape((3,nI*nJ*nK)).transpose(), order='C')

        ### old (wrong order) ###
        #axes = make_axes(xlims, ylims, zlims, gridspacing)
        #grid = make_grid(axes)

        gridspacings.append(gridspacing)
        allgrids.append(grid)

        counter += 1

    TreeGrid = np.array(allgrids)
    reconstructed_points = np.vstack(allgrids)
    assert(np.all( reconstructed_points.reshape(TreeGrid.shape) == TreeGrid ))
    print('TreeGrid.shape == ' + str(TreeGrid.shape))

    if testOUTvsRECON:
        print('individually reconstructed points from .tree file:')
        print(reconstructed_points)
        print('is reconstructed points and .out points the same in the same order:')
        print(np.all(points == reconstructed_points))
        print('is reconstructed points and .out points the same in potentially different order:')
        print(isReordering(points, reconstructed_points))

    ind = Reorder(points, reconstructed_points)
    assert(np.all( points[ind,:] == reconstructed_points ))
    assert(np.all( points[ind,:].reshape((len(blockused_), nI*nJ*nK, 3)) == TreeGrid ))

    if save:
        np.save(fname+'.out-indices_for_reordering.npy', ind)
    return ind



if __name__ == '__main__':
    ind = generate_indices("../../Batsrus.jl-master/3d__var_1_n00002500")
    print(ind)
