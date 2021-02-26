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


def get_blocks_used(dTREE):
    blockused_ = [iNode for iNode in range(1, dTREE['iTree_IA'].shape[1]+1) if dTREE['iTree_IA'][FortEqu(Status_), FortEqu(iNode)] == Used_]
    return blockused_


# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951
def get_tree_position(iNode, dTREE, returnall=False):
    # Calculate normalized position of the edges of node inode.
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





'''
# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975
def find_tree_node(CoordIn_D):
    ### todo: generalize later:
    nDimAmr = 3
    iRatio = 2
    jRatio = 2
    kRatio = 2
    ###

    ############ from SWMF/GM/BATSRUS/srcBATL/BATL_size_orig.f90 
    # Maximum dimensionality of grid is 3 (cannot be modified)
    MaxDim = 3

    # Indexes of AMR dimensions.
    # The magic formulas should be correct from 1 to nDimAmr.
    iDimAmrTmp_D = np.empty((MaxDim,))
    iDimAmrTmp_D = np.array([1 + (2-iRatio)*(3-jRatio), 6-iRatio-jRatio, 3 ], dtype=np.int8)
    print('iDimAmrTmp_D = ' + str(iDimAmrTmp_D))

    iDimAmr_D = iDimAmrTmp_D[0:nDimAmr]
    ############

    nLevelMin, nLevelMax = order_tree()

    # Find the node that contains a point. The point coordinates should
    # be given in generalized coordinates normalized to the domain size:
    # CoordIn_D = (CoordOrig_D - CoordMin_D)/DomainSize_D

    ##real, intent(in):: CoordIn_D(MaxDim)
    ##integer, intent(out):: iNode

    ##real :: Coord_D(MaxDim)
    ##integer :: iLevel, iChild
    ##integer :: iRoot_D(MaxDim), iCoord_D(nDimAmr), iBit_D(nDimAmr)

    # Scale coordinates so that 1 <= Coord_D <= nRoot_D+1
    #--------------------------------------------------------------------------
    Coord_D = 1.0 + nRoot_D*np.maximum(0.0, np.minimum(1.0, CoordIn_D))

    # Get root node index
    iRoot_D = np.minimum( np.array(Coord_D,dtype=np.int32),  nRoot_D )

    # Root node indexes are ordered
    iNode = iRoot_D[0] + nRoot_D[0]*((iRoot_D[1]-1) + nRoot_D[1]*(iRoot_D[2]-1))

    if(iTree_IA[FortEqu(Status_),FortEqu(iNode)] == Used_): 
        return iNode

    # Get normalized coordinates within root node and scale it up
    # to the largest resolution: 0 <= iCoord_D <= MaxCoord_I(nLevelMax)-1
    iCoord_D = np.minimum(MaxCoord_I[FortEqu(nLevelMax)] - 1,
                        np.array((Coord_D[FortEqu(iDimAmr_D)] - iRoot_D[FortEqu(iDimAmr_D)])*MaxCoord_I[FortEqu(nLevelMax)],
                                                        dtype=np.int32) 
                        )

    
    # Go down the tree using bit information
    for iLevel in range(nLevelMax-1,-1,-1):
        # Get the binary bits based on the coordinates
        iBit_D = ibits(iCoord_D, iLevel, 1)#???FortEqu
        # Construct child index as iChild = Sum Bit_i*2**i
        # The powers of 2 are stored in MaxCoord_I
        iChild = np.sum(iBit_D*MaxCoord_I[0:nDimAmr]) + Child1_
        iNode  = iTree_IA[FortEqu(iChild),FortEqu(iNode)]

        if iTree_IA[FortEqu(Status_),FortEqu(iNode)] == Used_:
            return iNode
    

    # Did not find the point so set iNode as unset
    iNode = Unset_
    return iNode
'''
