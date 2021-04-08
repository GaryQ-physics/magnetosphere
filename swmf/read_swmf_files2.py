import numpy as np
from swmf_constants import Used_,Status_,Level_,Parent_,Child0_,Child1_,Coord1_,CoordLast_,ROOTNODE_
import util
from named_var_indexes import nVarNeeded, index2str

def F2P(fortran_index):
    return fortran_index - 1

def P2F(python_index):
    return python_index + 1


def read_info_file(filetag):
    cache = {'filetag' : filetag}
    with open(filetag+'.info','r') as f:
        for line in f.readlines():
            if line == '\n' : continue
            if line[0] == '#': continue
            splt = line.split()
            if len(splt) == 2:
                cache[splt[1]] = splt[0]
    return cache

def read_tree_file(filetag):
    # first read info file
    cache = read_info_file(filetag)
    cache['nDim'] = int(cache['nDim'])
    cache['nI'] = int(cache['BlockSize1'])
    cache['nJ'] = int(cache['BlockSize2'])
    cache['nK'] = int(cache['BlockSize3'])

    ## Loading AMR tree
    # directly read bytes
    #f = open(filetag+".tree", 'rb')
    #filebytes = np.fromfile(f, dtype=np.int32)
    #f.close()

    try:
        # use scipy FortranFile
        from scipy.io import FortranFile
        ff = FortranFile(filetag+".tree", 'r')
        if True:
            nDim, nInfo, nNode = ff.read_reals(dtype=np.int32)
            iRatio_D = ff.read_reals(dtype=np.int32) # Array of refinement ratios
            nRoot_D = ff.read_reals(dtype=np.int32) # The number of root nodes in all dimension
            iTree_IA = ff.read_reals(dtype=np.int32).reshape((nInfo,nNode), order='F')
        else:
            nDim, nInfo, nNode = ff.read_ints(dtype='i4')
            iRatio_D = ff.read_ints(dtype='i4') # Array of refinement ratios
            nRoot_D = ff.read_ints(dtype='i4') # The number of root nodes in all dimension
            iTree_IA = ff.read_ints(dtype='i4').reshape((nInfo,nNode), order='fortran')
    except:
        raise RuntimeWarning ("scipy.io.FortranFile didnt work")
        # use fortranfile
        from fortranfile import FortranFile
        ff = FortranFile(filetag+".tree") # read or write ???
        nDim, nInfo, nNode = ff.readInts()
        iRatio_D = ff.readInts() # Array of refinement ratios
        nRoot_D = ff.readInts() # The number of root nodes in all dimension
        iTree_IA = ff.readInts().reshape((nInfo,nNode), order='fortran')

    ########################### check_thing_work #######################
    assert(cache['nDim'] == nDim)
    assert(iTree_IA.shape[1] == nNode)
    # Maximum number of ghost cells set by Config.pl script.
    # Valid values are 0,1,2,3,4,5
    nG = 2
    # Refinement ratios in the 3 dimensions. Either 1 or 2.
    # The values are set by the Config.pl script.
    iRatio, jRatio, kRatio = min(2, cache['nI']), min(2, cache['nJ']), min(2, cache['nK'])
    # Number of dimensions in which grid adaptation is done
    nDimAmr = iRatio + jRatio + kRatio - 3
    assert(nDimAmr == nDim)
    assert(nDim == 3)
    # Number of children per node
    nChild = 2**nDimAmr
    assert(np.isfortran(iTree_IA))
    assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))
    ####################################################################

    cache['nInfo']    = nInfo
    cache['nNode']    = nNode
    cache['iRatio_D'] = iRatio_D
    cache['nRoot_D']  = nRoot_D
    cache['iTree_IA'] = iTree_IA

    cache['xGlobalMin'] = float(cache['Coord1Min'])
    cache['yGlobalMin'] = float(cache['Coord2Min'])
    cache['zGlobalMin'] = float(cache['Coord3Min'])
    cache['xGlobalMax'] = float(cache['Coord1Max'])
    cache['yGlobalMax'] = float(cache['Coord2Max'])
    cache['zGlobalMax'] = float(cache['Coord3Max'])

    return cache


def read_out_file(filetag):
    pass


# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951 with substitutions
def get_tree_position(iNode, cache, returnall=False):
    '''
    Calculate normalized position of the edges of node iNode.
    Zero is at the minimum boundary of the grid, one is at the max boundary
    '''
    iTree_IA = cache['iTree_IA']
    nRoot_D  = cache['nRoot_D']
    iRatio_D = cache['iRatio_D']

    iLevel = iTree_IA[F2P(Level_), F2P(iNode)]

    MaxIndex_D = ((2**(iLevel)-1)*(iRatio_D-1) + 1)*nRoot_D
    # Note: in the common case of iRatio_D=[2,2,2] and nRoot_D=[1,1,1]:
    # MaxIndex_D[all] = ((2**(iLevel)-1)*(2-1) + 1)*1
    #                 = 2**iLevel

    assert(MaxIndex_D.shape == (3,))
    assert(np.all( MaxIndex_D == 2**iLevel ))
    # note that if gridspacing = (256./8.)*0.5**iLevel, then gridspacing*MaxIndex_D[all] == 256./8. == 32. )

    # Convert to real by adding -1.0 or 0.0 for the two edges, respectively
    block_coords = iTree_IA[F2P(Coord1_):F2P(CoordLast_)+1,F2P(iNode)] # not seperatly defined in BATL_tree.f90
    PositionMin_D = (block_coords - 1.0)/MaxIndex_D
    PositionMax_D = (block_coords + 0.0)/MaxIndex_D

    if returnall:
        return PositionMin_D, PositionMax_D, MaxIndex_D, BlockCoord_D
    else:
        return PositionMin_D, PositionMax_D # what was returned in original


def get_physical_dimensions(iNode, cache, returnCenters=False):
    x_start = cache['xGlobalMin']
    y_start = cache['yGlobalMin']
    z_start = cache['zGlobalMin']
    x_range = cache['xGlobalMax'] - cache['xGlobalMin']
    y_range = cache['yGlobalMax'] - cache['yGlobalMin']
    z_range = cache['zGlobalMax'] - cache['zGlobalMin']

    iLevel = cache['iTree_IA'][F2P(Level_), F2P(iNode)]
    assert(cache['nI'] == cache['nJ'] == cache['nK'])
    assert(x_range == y_range == z_range)
    gridspacing = (x_range/cache['nI'])*0.5**iLevel

    PositionMin_D, PositionMax_D = get_tree_position(iNode, cache)
    xmin = x_range*(PositionMin_D[0]) + x_start
    ymin = y_range*(PositionMin_D[1]) + y_start
    zmin = z_range*(PositionMin_D[2]) + z_start
    xmax = x_range*(PositionMax_D[0]) + x_start
    ymax = y_range*(PositionMax_D[1]) + y_start
    zmax = z_range*(PositionMax_D[2]) + z_start

    xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
    ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
    zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)

    xminmax = (xmin, xmax)
    yminmax = (ymin, ymax)
    zminmax = (zmin, zmax)

    if returnCenters:
        return xlims, ylims, zlims, gridspacing
    else:
        return xminmax, yminmax, zminmax, gridspacing


# supposed to reproduce SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975, but differently
def find_tree_node(point, cache):
    iTree_IA = cache['iTree_IA']

    xin = cache['xGlobalMin'] <= point[0] <= cache['xGlobalMax']
    yin = cache['yGlobalMin'] <= point[1] <= cache['yGlobalMax']
    zin = cache['zGlobalMin'] <= point[2] <= cache['zGlobalMax']

    if not (xin and yin and zin): 
        raise RuntimeError ('point out of simulation volume')

    iNode = ROOTNODE_
    while(True):
        if Used_ == iTree_IA[F2P(Status_), F2P(iNode)]:
            break

        for j in range(8):
            child = iTree_IA[F2P(Child1_+j), F2P(iNode)]
            xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(child, cache, returnCenters=False)

            xin = xminmax[0] <= point[0] <= xminmax[1]
            yin = yminmax[0] <= point[1] <= yminmax[1]
            zin = zminmax[0] <= point[2] <= zminmax[1]

            if xin and yin and zin: 
                iNode = child
                break

    return iNode


def read_all(filetag):
    def Reorder(arr1, arr2):
        ''' arr1 and arr2 are (N,3) arrays of for N  3d-points 
            if arr1 and arr2 contain the same points, but in different
            order, returns array of indices 'ind' such that
            arr1[ind,:] == arr2
        '''
        arr1 = np.array(arr1); arr2 = np.array(arr2)
        assert(len(arr1.shape) == len(arr2.shape) == 2)
        if arr1.shape != arr2.shape: return False

        sort1 = np.lexsort((arr1[:,0],arr1[:,1],arr1[:,2]))  
        sort2 = np.lexsort((arr2[:,0],arr2[:,1],arr2[:,2]))
        #undo1 = np.argsort(sort1)
        undo2 = np.argsort(sort2)

        if not np.all(arr2 == (arr1[sort1[undo2],:])):
            raise RuntimeError ('arrays inputed arent reordering of each other')
        return sort1[undo2]

    if True:
        import spacepy.pybats.bats as bats
        import read_swmf_files as rswmf
        data = bats.Bats2d(filetag + ".out")
        header = "R R R Mp/cc km/s km/s km/s J/m3 nT nT nT nT nT nT nPa uA/m2 uA/m2 uA/m2 --"
        assert(header == data.meta['header'].strip())
    else:
        data = read_out_file(filetag) # TODO : eventually
    cache = read_tree_file(filetag)

    # in what follows:
    #  the P in iNodeP and iBlockP stands for python like indexing (as oposed to fortran)
    #  
    #  iNodeP indexes all nodes of the tree, from 0 to nNode-1,
    #  and thus the "iNode" to be used in the other functions is simply iNodeP+1, or P2F(iNodeP)
    # 
    #  iBlockP indexes all the blocks used, from 0 to nBlock-1.
    #  There is one for each node with a status of used. 
    #  Note, nBlock*nI*nJ*nK = total number of batsrus cells (npts)
    points = np.column_stack([data['x'],data['y'],data['z']])
    npts = points.shape[0]
    nI, nJ, nK = cache['nI'], cache['nJ'], cache['nK']
    nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1
    nNode = cache['nNode']

    iBlockP = 0
    block2node = -np.ones((nBlock,), dtype=int)#!!
    node2block = -np.ones((nNode,), dtype=int)#!!
    reconstructed_points = np.nan*np.empty((npts,3), dtype=np.float32)
    for iNodeP in range(nNode):
        if cache['iTree_IA'][F2P(Status_), iNodeP] == Used_:
            block2node[iBlockP] = iNodeP
            node2block[iNodeP] = iBlockP
            xlims, ylims, zlims, gridspacing = rswmf.get_physical_dimensions(P2F(iNodeP), cache, returnCenters=True)
            assert(xlims[1]-xlims[0] == (nI-1)*gridspacing)
            assert(ylims[1]-ylims[0] == (nJ-1)*gridspacing)
            assert(zlims[1]-zlims[0] == (nK-1)*gridspacing)

            grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                            ylims[0]:ylims[1]+gridspacing:gridspacing,
                            zlims[0]:zlims[1]+gridspacing:gridspacing ]
            grid = np.array(grid.reshape((3,nI*nJ*nK)).transpose(), order='C')

            #### equivalent too ####
            #grid = np.empty((nI,nJ,nK,3))
            #for i in range(nI):
            #    for j in range(nJ):
            #        for k in range(nK):
            #            grid[i,j,k, 0] = xlims[0]+gridspacing*i
            #            grid[i,j,k, 1] = ylims[0]+gridspacing*j
            #            grid[i,j,k, 2] = zlims[0]+gridspacing*k
            #grid = grid.reshape((nI*nJ*nK,3))
            ############

            start = iBlockP*nI*nJ*nK
            end = (iBlockP+1)*nI*nJ*nK
            reconstructed_points[start:end,:] = grid

            iBlockP += 1

    ind = Reorder(points, reconstructed_points)

    DataArray = np.empty((nVarNeeded, nBlock, nI, nJ, nK), dtype=np.float32); DataArray[:,:,:,:,:] = np.nan
    for index in range(nVarNeeded):
        DataArray[index,:,:,:,:] = data[index2str[index]][ind].reshape((nBlock, nI, nJ, nK))
    del data, ind

    cache['DataArray'] = DataArray
    cache['block2node'] = block2node
    cache['node2block'] = node2block
    cache['nBlock'] = nBlock
    return cache


def find_index(point):
    pass


def interpolate(filetag, point, var='p', cache=None, debug=False):
    """
    arguments:
        filetag:
            string with the name of the swmf output files,
            including the path, not includint the extension.
            e.g. '/home/user/data/3d__var_3_e20010101-010000-000'
        point(s):

        var (optional):
            string for the swmf variable name. Default: var='p'
    returns:
    """
    #TODO: vectorize find_tree_node for points (N,3) (possibly with jitFORTRAN).
    #      then get rid of get_physical_dimensions() call in this function
    #      gridpacing and minmax can be found direcly from block_data.
    #      then maybe this function can be vectorized as well for points (N,3)
    # this function should maybe go in different file also

    if cache is None:
        cache = read_all(filetag)
    else:
        assert(cache['filetag'] == filetag)

    def getvar(_var, iNode, i,j,k):
        return cache['DataArray'][_var,:,:,:,:][node2block[F2P(iNode)],i,j,k]

    iNode = find_tree_node(point, cache)

    # get the gridspacing in x,y,z
    gridspacingX = getvar('x', iNode, 1,0,0) - getvar('x', iNode, 0,0,0)
    gridspacingY = getvar('y', iNode, 0,1,0) - getvar('y', iNode, 0,0,0)
    gridspacingZ = getvar('z', iNode, 0,0,1) - getvar('z', iNode, 0,0,0)

    # i0 is s.t. the highest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
    i0 = int(np.floor( (point[0] - getvar('x', iNode, 0, 0, 0))/gridspacingX ))
    j0 = int(np.floor( (point[1] - getvar('y', iNode, 0, 0, 0))/gridspacingY ))
    k0 = int(np.floor( (point[2] - getvar('z', iNode, 0, 0, 0))/gridspacingZ ))
    if debug: print(getvar('x', iNode, 0, 0, 0),getvar('y', iNode, 0, 0, 0),getvar('z', iNode, 0, 0, 0))
    if debug: print(i0,j0,k0)
    #i1 = i0+1 is the lowest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i1,:,:]  is still greater than point[0]

    # together, i0 and i1 form the upper and lower bounds for a linear interpolation in x
    # likewise for j0,j1,y  and k0,k1,z

    #TODO implement better interpolation at ends of block.
    # This method effectively makes it nearest neighbor at the ends
    if i0 == -1:
        i0 = 0
        i1 = 0
        if debug: print('edge case')
    elif i0 == cache['nI']-1:
        i1 = cache['nI']-1
        if debug: print('edge case')
    else:
        i1 = i0 + 1

    if j0 == -1:
        j0 = 0
        j1 = 0
        if debug: print('edge case')
    elif j0 == cache['nJ']-1:
        j1 = cache['nJ']-1
        if debug: print('edge case')
    else:
        j1 = j0 + 1

    if k0 == -1:
        k0 = 0
        k1 = 0
        if debug: print('edge case')
    elif k0 == cache['nK']-1:
        k1 = cache['nK']-1
        if debug: print('edge case')
    else:
        k1 = k0 + 1

    # all together i0,i1,j0,ect... form a cube of side length "gridpacing" to do trililear interpolation within
    # define xd as the distance along x of point within that cube, in units of "gridspacing"
    xd = (point[0] - getvar('x', iNode, i0, 0 , 0 ) )/gridspacingX
    yd = (point[1] - getvar('y', iNode, 0 , j0, 0 ) )/gridspacingY
    zd = (point[2] - getvar('z', iNode, 0 , 0 , k0) )/gridspacingZ

    if debug: print(xd,yd,zd)
    if debug: print(getvar('x', iNode,  i0, j0, k0))
    if debug: print(getvar('y', iNode,  i0, j0, k0))
    if debug: print(getvar('z', iNode,  i0, j0, k0))
    if debug: print('hellothere')
    if debug: print((iNode, node2block[F2P(iNode)], i0, j0, k0))


    #https://en.wikipedia.org/wiki/Trilinear_interpolation
    c000 = getvar(var, iNode,  i0, j0, k0)
    c001 = getvar(var, iNode,  i0, j0, k1)
    c010 = getvar(var, iNode,  i0, j1, k0)
    c100 = getvar(var, iNode,  i1, j0, k0)
    c011 = getvar(var, iNode,  i0, j1, k1)
    c110 = getvar(var, iNode,  i1, j1, k0)
    c101 = getvar(var, iNode,  i1, j0, k1)
    c111 = getvar(var, iNode,  i1, j1, k1)
    if debug: print(c000)
    if debug: print(c001,c010,c100,c100,c011,c110,c101,c111)

    c00 = c000*(1.-xd) + c100*xd
    c01 = c001*(1.-xd) + c101*xd
    c10 = c010*(1.-xd) + c110*xd
    c11 = c011*(1.-xd) + c111*xd

    c0 = c00*(1.-yd) + c10*yd
    c1 = c01*(1.-yd) + c11*yd

    c = c0*(1.-zd) + c1*zd
    if debug: print(c)
    return c



if __name__ == '__main__':
    point = np.array([-220., -124.,  -92.])
    print(interpolate("/home/gary/temp/3d__var_3_e20031120-070000-000",point,var='p', debug=True))

    point = np.array([1., 0.,  0.])
