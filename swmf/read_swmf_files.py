import numpy as np

def FortEqu(fortranindex):
    return fortranindex - 1

# Possible values for the status variable
Unset_     = -100 # index for unset values (that are otherwise larger)
Unused_      = -1 # unused block (not a leaf)
Refine_      = -2 # parent block to be refined
DontCoarsen  = -3 # block not to be coarsened
Coarsen_     = -4 # child block to be coarsened
Used_        =  1 # currently used block (leaf)
RefineNew_   =  2 # child block to be refined
Refined_     =  3 # refined child block
CoarsenNew_  =  4 # parent block to be coarsened
Coarsened_   =  5 # coarsened parent block

# Deepest AMR level relative to root nodes (limited by 32 bit integers)
MaxLevel = 30

'''from SWMF/GM/BATSRUS/srcBATL_/BATL_tree.f90
  ! The maximum integer coordinate for a given level below root nodes
  ! Implied do loop was not understooed by the pgf90 compiler, so list them
  integer, parameter, public :: MaxCoord_I(0:MaxLevel) = &
       [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, &
       16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, &
       4194304, 8388608, 16777216, 33554432, 67108864, 134217728, &
       268435456, 536870912, 1073741824 ]
'''#why are they indexing this one from 0?

MaxCoord_I =  [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 
                4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
                268435456, 536870912, 1073741824 ]

# check the copied hardcoded things are what I think they are
assert(len(MaxCoord_I) == MaxLevel+1)
assert(MaxCoord_I == [2**i for i in range(len(MaxCoord_I))])

# Named indexes of iTree_IA
Status_   =  1
Level_    =  2 # grid level
Proc_     =  3 # processor index
Block_    =  4 # block index
MinLevel_ =  5 # minimum level allowed
MaxLevel_ =  6 # maximum level allowed
Coord0_   =  6 # equal to Coord1_-1
Coord1_   =  7 # coordinate of node in 1st dimension
Coord2_   =  8 # coordinate of node in 2nd dimension
Coord3_   =  9 # coordinate of node in 3rd dimension
CoordLast_=  9 # Coord0_ + MaxDim (?)
Parent_   = 10 # Parent_ must be equal to Child0_
Child0_   = 10 #
Child1_   = Child0_ + 1
#ChildLast_= Child0_ + nChild

'''
the status of the node (used, unused, to be refined, to be coarsened, etc.);
the current, the maximum allowed and minimum allowed AMR levels for this node;
the three integer coordinates with respect to the whole grid;
the index of the parent node (if any);
the indexes of the children nodes (if any);
the processor index where the block is stored for active nodes;
the local block index for active nodes.
'''

# my added named constants
iRootNode = 1

def read_info_file(filetag):
    dINFO={}
    with open(filetag+'.info','r') as f:
        for line in f.readlines():
            if line == '\n' : continue
            if line[0] == '#': continue
            splt = line.split()
            if len(splt) == 2:
                dINFO[splt[1]] = splt[0]
    return dINFO


def read_tree_file(filetag):
    # first read info file
    dINFO=read_info_file(filetag)
    nDim = int(dINFO['nDim'])
    nI = int(dINFO['BlockSize1'])
    nJ = int(dINFO['BlockSize2'])
    nK = int(dINFO['BlockSize3'])

    ## Loading AMR tree
    # directly read bytes
    f = open(filetag+".tree", 'rb')
    filebytes = np.fromfile(f, dtype=np.int32)
    f.close()

    usesp = True
    if usesp:
        # use scipy FortranFile
        from scipy.io import FortranFile
        ff = FortranFile(filetag+".tree", 'r')
        if True:
            nDim, nInfo, nNode = ff.read_reals(dtype=np.int32)
            iRatio_D = ff.read_reals(dtype=np.int32) # Array of refinement ratios
            nRoot_D = ff.read_reals(dtype=np.int32) # The number of root nodes in all dimension
            iTree_IA = ff.read_reals(dtype=np.int32).reshape((nInfo,nNode), order='F')
            #iTree_IA = ff.read_reals(dtype=np.int32).reshape((nNode,nInfo)).transpose()
        else:
            nDim, nInfo, nNode = ff.read_ints(dtype='i4')
            iRatio_D = ff.read_ints(dtype='i4') # Array of refinement ratios
            nRoot_D = ff.read_ints(dtype='i4') # The number of root nodes in all dimension
            iTree_IA = ff.read_ints(dtype='i4').reshape((nInfo,nNode), order='fortran')
    else:
        # use fortranfile
        from fortranfile import FortranFile
        ff = FortranFile(filetag+".tree") # read or write ???
        nDim, nInfo, nNode = ff.readInts()
        iRatio_D = ff.readInts() # Array of refinement ratios
        nRoot_D = ff.readInts() # The number of root nodes in all dimension
        iTree_IA = ff.readInts().reshape((nInfo,nNode), order='fortran')

    ########################### check_thing_work #######################
    # Maximum number of ghost cells set by Config.pl script.
    # Valid values are 0,1,2,3,4,5
    nG = 2
    # Refinement ratios in the 3 dimensions. Either 1 or 2.
    # The values are set by the Config.pl script.
    iRatio, jRatio, kRatio = min(2, nI), min(2, nJ), min(2, nK)
    # Number of dimensions in which grid adaptation is done
    nDimAmr = iRatio + jRatio + kRatio - 3
    assert(nDimAmr == nDim)
    assert(nDim == 3)
    # Number of children per node
    nChild = 2**nDimAmr
    assert(np.isfortran(iTree_IA))
    assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))
    ####################################################################

    dTREE = {'iTree_IA': iTree_IA, 'nRoot_D': nRoot_D, 'iRatio_D': iRatio_D,
             'nDims': 3, 'nI': nI, 'nJ': nJ, 'nK': nK }
    dTREE['xGlobalMin'] = float(dINFO['Coord1Min'])
    dTREE['yGlobalMin'] = float(dINFO['Coord2Min'])
    dTREE['zGlobalMin'] = float(dINFO['Coord3Min'])
    dTREE['xGlobalMax'] = float(dINFO['Coord1Max'])
    dTREE['yGlobalMax'] = float(dINFO['Coord2Max'])
    dTREE['zGlobalMax'] = float(dINFO['Coord3Max'])

    return dTREE


def read_out_file(filetag):
    pass


def read_all(filetag):
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

    if True:
        import spacepy.pybats.bats as bats
        import read_swmf_files as rswmf
        data = bats.Bats2d(filetag + ".out")
    else:
        data = read_out_file(filetag) # eventually
    dTREE = read_tree_file(filetag)

    points = np.column_stack([data['x'],data['y'],data['z']])
    npts = points.shape[0]
    nI, nJ, nK = dTREE['nI'], dTREE['nJ'], dTREE['nK']
    nBlocks = int(npts/(nI*nJ*nK)) ##TODO: Check if proper terminology
    nNode = dTREE['iTree_IA'].shape[1] ##TODO: Check if proper terminology

    counter = 0
    blockused_ = -np.ones((nBlocks,), dtype=int)
    reverse_blockused_ = -np.ones((nNode,), dtype=int)
    reconstructed_points = np.nan*np.empty((npts,3))
    for iNode in range(1, nNode+1):
        if dTREE['iTree_IA'][FortEqu(Status_), FortEqu(iNode)] == Used_:
            blockused_[counter] = iNode
            reverse_blockused_[FortEqu(iNode)] = counter 
            xlims, ylims, zlims, gridspacing = rswmf.get_physical_dimensions(iNode, dTREE, returnCenters=True)
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

            start = counter*nI*nJ*nK
            end = (counter+1)*nI*nJ*nK
            reconstructed_points[start:end,:] = grid

            counter += 1

    ind = Reorder(points, reconstructed_points)
    #assert(np.all( points[ind,:] == reconstructed_points ))
    #np.save(fname+'.out-indices_for_reordering.npy', ind)
    #print(npts)
    #print(int(npts/(nI*nJ*nK)))
    #print(dTREE['iTree_IA'].shape)
    #print(len(blockused_)) #11516
    #print(len(reverse_blockused_)) #13161
    #print(reverse_blockused_)
    #print(blockused_)
    #print(ind)
    return [data, dTREE, ind, blockused_, reverse_blockused_]

# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951
def get_tree_position(iNode, dTREE, returnall=False):
    # Calculate normalized position of the edges of node iNode.
    # Zero is at the minimum boundary of the grid, one is at the max boundary
    iTree_IA = dTREE['iTree_IA']
    nRoot_D = dTREE['nRoot_D']
    iRatio_D = dTREE['iRatio_D']

    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]

    # For non-AMR directions MaxIndex_D = nRoot_D
    # For AMR     directions MaxIndex_D = nRoot_D*MaxCoord_I(iLevel)
    MaxIndex_D = ((MaxCoord_I[iLevel]-1)*(iRatio_D-1) + 1)*nRoot_D # note, MaxCoord_I is indexed by 0 even in fortran code (so no need for FortEqu)
    # in this case:
    # MaxIndex_D[all] = ((2**(iLevel)-1)*(2-1) + 1)*1
    #                 = 2**iLevel
    assert(MaxIndex_D.shape == (3,))
    assert(np.all( MaxIndex_D == 2**iLevel ))
    # note that if gridspacing = (256./8.)*0.5**iLevel, then gridspacing*MaxIndex_D[all] == 256./8. == 32. )

    # Convert to real by adding -1.0 or 0.0 for the two edges, respectively
    BlockCoord_D = iTree_IA[FortEqu(Coord1_):FortEqu(CoordLast_)+1,FortEqu(iNode)] # not seperated defined in original. rename???
    PositionMin_D = (BlockCoord_D - 1.0)/MaxIndex_D
    PositionMax_D = (BlockCoord_D + 0.0)/MaxIndex_D

    if returnall:
        return PositionMin_D, PositionMax_D, MaxIndex_D, BlockCoord_D
    else:
        return PositionMin_D, PositionMax_D # what was returned in original


def get_physical_dimensions(iNode, dTREE, returnCenters=False):
    x_start = dTREE['xGlobalMin']
    y_start = dTREE['yGlobalMin']
    z_start = dTREE['zGlobalMin']
    x_range = dTREE['xGlobalMax'] - dTREE['xGlobalMin']
    y_range = dTREE['yGlobalMax'] - dTREE['yGlobalMin']
    z_range = dTREE['zGlobalMax'] - dTREE['zGlobalMin']

    iLevel = dTREE['iTree_IA'][FortEqu(Level_), FortEqu(iNode)]
    assert(dTREE['nI'] == dTREE['nJ'] == dTREE['nK'])
    assert(x_range == y_range == z_range)
    gridspacing = (x_range/dTREE['nI'])*0.5**iLevel

    PositionMin_D, PositionMax_D = get_tree_position(iNode, dTREE)
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
def find_tree_node(point, dTREE):
    iTree_IA = dTREE['iTree_IA']

    xin = dTREE['xGlobalMin'] <= point[0] <= dTREE['xGlobalMax']
    yin = dTREE['yGlobalMin'] <= point[1] <= dTREE['yGlobalMax']
    zin = dTREE['zGlobalMin'] <= point[2] <= dTREE['zGlobalMax']

    if not (xin and yin and zin): 
        raise RuntimeError ('point out of simulation volume')

    iNode = iRootNode
    while(True):
        if Used_ == iTree_IA[FortEqu(Status_), FortEqu(iNode)]:
            break

        for j in range(8):
            child = iTree_IA[FortEqu(Child1_+j), FortEqu(iNode)]
            xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(child, dTREE, returnCenters=False)

            xin = xminmax[0] <= point[0] <= xminmax[1]
            yin = yminmax[0] <= point[1] <= yminmax[1]
            zin = zminmax[0] <= point[2] <= zminmax[1]

            if xin and yin and zin: 
                iNode = child
                break

    return iNode


def test(dTREE):
    assert(dTREE['iTree_IA'][FortEqu(Level_), FortEqu(iRootNode)] == 0)
    xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(iRootNode, dTREE, returnCenters=False)
    assert(xminmax == (dTREE['xGlobalMin'], dTREE['xGlobalMax']))
    assert(yminmax == (dTREE['yGlobalMin'], dTREE['yGlobalMax']))
    assert(zminmax == (dTREE['zGlobalMin'], dTREE['zGlobalMax']))


def interpolate(filetag, point, var='p', debug=False):
    #TODO: vectorize find_tree_node for points (N,3) (possibly with jitFORTRAN).
    #      then get rid of get_physical_dimensions() call in this function
    #      gridpacing and minmax can be found direcly from block_data.
    #      then maybe this function can be vectorized as well for points (N,3)
    # this function should maybe go in different file also

    data, dTREE, ind, blockused_, reverse_blockused_ = read_all(filetag)
    nI, nJ, nK = dTREE['nI'], dTREE['nJ'], dTREE['nK']

    def getvar(_var, iNode, i,j,k):
        return data[_var][ind].reshape(len(blockused_),nI,nJ,nK)[reverse_blockused_[FortEqu(iNode)],i,j,k]

    iNode = find_tree_node(point, dTREE)
    xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(iNode, dTREE, returnCenters=False)
    #if debug: print(xminmax, yminmax, zminmax, gridspacing)
    #if debug: print(getvar('x',iNode,0,0,0))
    #if debug: print(getvar('y',iNode,0,0,0))
    #if debug: print(getvar('z',iNode,0,0,0))
    #if debug: print(getvar('x',iNode,1,0,0))
    #if debug: print(getvar('y',iNode,1,0,0))
    #if debug: print(getvar('z',iNode,1,0,0))
    #if debug: print(getvar('x',iNode,0,1,0))
    #if debug: print(getvar('y',iNode,0,1,0))
    #if debug: print(getvar('z',iNode,0,1,0))
    #if debug: print(getvar('x',iNode,0,0,2))
    #if debug: print(getvar('y',iNode,0,0,2))
    #if debug: print(getvar('z',iNode,0,0,2))

    # i0 is s.t. the highest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
    i0 = int(np.floor( (point[0] - getvar('x', iNode, 0, 0, 0))/gridspacing ))
    j0 = int(np.floor( (point[1] - getvar('y', iNode, 0, 0, 0))/gridspacing ))
    k0 = int(np.floor( (point[2] - getvar('z', iNode, 0, 0, 0))/gridspacing ))

    #i1 = i0+1 is the lowest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i1,:,:]  is still greater than point[0]

    # together, i0 and i1 form the upper and lower bounds for a linear interpolation in x
    # likewise for j0,j1,y  and k0,k1,z

    #TODO implement better interpolation at ends of block.
    # This effectively makes it nearest neighbor at the ends
    if i0 == -1:
        i0 = 0
        i1 = 0
        if debug: print('edge case')
    elif i0 == dTREE['nI']-1:
        i1 = dTREE['nI']-1
        if debug: print('edge case')
    else:
        i1 = i0 + 1

    if j0 == -1:
        j0 = 0
        j1 = 0
        if debug: print('edge case')
    elif j0 == dTREE['nJ']-1:
        j1 = dTREE['nJ']-1
        if debug: print('edge case')
    else:
        j1 = j0 + 1

    if k0 == -1:
        k0 = 0
        k1 = 0
        if debug: print('edge case')
    elif k0 == dTREE['nK']-1:
        k1 = dTREE['nK']-1
        if debug: print('edge case')
    else:
        k1 = k0 + 1

    # all together i0,i1,j0,ect... form a cube of side length "gridpacing" to do trililear interpolation within
    # define xd as the distance along x of point within that cube, in units of "gridspacing"
    xd = (point[0] - getvar('x', iNode, i0, 0 , 0 ) )/gridspacing
    yd = (point[1] - getvar('y', iNode, 0 , j0, 0 ) )/gridspacing
    zd = (point[2] - getvar('z', iNode, 0 , 0 , k0) )/gridspacing

    if debug: print(xd,yd,zd)
    if debug: print(getvar('x', iNode,  i0, j0, k0))
    if debug: print(getvar('y', iNode,  i0, j0, k0))
    if debug: print(getvar('z', iNode,  i0, j0, k0))

    #TODO make better?
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
    print(interpolate("/home/gary/temp/3d__var_3_e20031120-070000-000",point,var='p'))
