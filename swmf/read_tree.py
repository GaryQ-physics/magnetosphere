import numpy as np
import struct

debug = False

def FortEqu(fortranindex):
    return fortranindex - 1

##
# Number of cells per block in each direction.
# These values are set by the Config.pl script.
# Set 1 for ignored directions!
nI, nJ, nK = 8, 8, 8

# Maximum number of ghost cells set by Config.pl script.
# Valid values are 0,1,2,3,4,5
nG = 2

# Refinement ratios in the 3 dimensions. Either 1 or 2.
# The values are set by the Config.pl script.
iRatio, jRatio, kRatio = min(2, nI), min(2, nJ), min(2, nK)

# Number of dimensions in which grid adaptation is done
nDimAmr = iRatio + jRatio + kRatio - 3

assert(nDimAmr == 3)

# Number of children per node
nChild = 2**nDimAmr

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
ChildLast_= Child0_ + nChild

'''
the status of the node (used, unused, to be refined, to be coarsened, etc.);
the current, the maximum allowed and minimum allowed AMR levels for this node;
the three integer coordinates with respect to the whole grid;
the index of the parent node (if any);
the indexes of the children nodes (if any);
the processor index where the block is stored for active nodes;
the local block index for active nodes.
'''


filetag = "/home/gary/temp/3d__var_3_e20031120-070000-000"



## Loading AMR tree
usesp = True

# directly read bytes
f = open(filetag+".tree", 'rb')
filebytes = np.fromfile(f, dtype=np.int32)
f.close()
if debug: print(filebytes[:20])

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


assert(np.isfortran(iTree_IA))
assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))

# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line358
nRoot = np.product(nRoot_D)


'''WRONG
headerInt = filebytes[0]
nDim = filebytes[1]
nInfo = filebytes[2]
nNode = filebytes[3]
iRatio_D = filebytes[4:4+nDim] # Array of refinement ratios
nRoot_D = filebytes[4+nDim:4+2*nDim]  # The number of root nodes in all dimension
iTree_IA = filebytes[4+2*nDim:4+2*nDim+nInfo*nNode].reshape((nInfo,nNode))
'''
if debug:
    print('nDim = ' + str(nDim))
    print('nInfo = ' + str(nInfo))
    print('nNode = ' + str(nNode))
    print('iRatio_D = ' + str(iRatio_D))
    print('nRoot_D = ' + str(nRoot_D))

    print(nDim)
    print(nInfo)
    print(nNode)
    print(iRatio_D)
    print(nRoot_D)
    print(iTree_IA.shape)

    print(iTree_IA)
    print(iTree_IA[:,0])
    print(iTree_IA[0,:])


# Get all the used blocks
allblocks_ = [iTree_IA[FortEqu(Status_), FortEqu(iNode)] for iNode in range(1, nNode+1)]
blockused_ = [iNode for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
nNodeUsed = len(blockused_) # I'm assuming this is the same as one in BATL_tree.f90

MinLevels = [iTree_IA[FortEqu(MinLevel_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
MaxLevels = [iTree_IA[FortEqu(MaxLevel_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Blocks = [iTree_IA[FortEqu(Block_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Levels = [iTree_IA[FortEqu(Level_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Procs = [iTree_IA[FortEqu(Proc_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]

if debug:
    print((256./8.)*0.5**max(Levels))
    print((256./8.)*0.5**min(Levels))
    print(len(blockused_)*8**3)


    print('#############\n\n')
    print(set(MinLevels))
    print(set(MaxLevels))
    print(set(Procs))
    print(max(Procs))
    print('#############\n\n')

#print('#############\n\n')
#print(set(allblocks))
#print(set(TEST)==set(range(837))) for Proc_
#print(blockused_)

# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951
def get_tree_position(iNode, returnall=False):
    # Calculate normalized position of the edges of node inode.
    # Zero is at the minimum boundary of the grid, one is at the max boundary

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

'''
  subroutine order_tree

    ! Set iNodeMorton_I and iMortonNode_A indirect index arrays according to
    ! 1. root node order
    ! 2. Morton ordering for each root node

    integer :: iNode, iRoot, jRoot, kRoot, iLevel
    !--------------------------------------------------------------------------
    nNode = nRoot
    iNode = 0
    iMorton = 0
    iNodeMorton_I(1:nNodeUsed) = Unset_
    iMortonNode_A(1:nNodeUsed) = Unset_
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)
             ! Root nodes are the first ones
             iNode = iNode + 1

             ! All root nodes are handled as if they were first child
             call order_children(iNode)
          end do
       end do
    end do

    nNodeUsed = iMorton

    ! Set min and max refinement levels
    nLevelMin = MaxLevel
    nLevelMax = 0
    do iMorton = 1, nNodeUsed
       iNode = iNodeMorton_I(iMorton)
       iLevel = iTree_IA(Level_,iNode)
       nLevelMin = min(iLevel, nLevelMin)
       nLevelMax = max(iLevel, nLevelMax)
    end do

  end subroutine order_tree
  !============================================================================
  recursive subroutine order_children(iNode)

    ! Recursively apply Morton ordering for nodes below a root block.
    ! Store result into iNodeMorton_I and iMortonNode_A using the global
    ! iMorton index.

    integer, intent(in) :: iNode
    integer :: iChild

    !--------------------------------------------------------------------------
    nNode = max(nNode, iNode)

    if(iTree_IA(Status_, iNode) >= Used_)then
       iMorton = iMorton + 1
       iNodeMorton_I(iMorton) = iNode
       iMortonNode_A(iNode)   = iMorton
    else
       do iChild = Child1_, ChildLast_
          call order_children(iTree_IA(iChild, iNode))
       end do
    end if

  end subroutine order_children
'''

#from https://bytes.com/topic/python/answers/158954-bitwise-operations-python
def ibits(i,pos,length):
    if isinstance(i, np.ndarray):
        i = i.ravel()
        ret = np.empty(i.shape, dtype=np.int32)
        for k in range(i.size):
            ret[k] = ibits(i[k],pos,length)
        return ret

    return (i >> pos) & ~(-1 << length)


def order_tree(): # from subroutine order_tree in BATL_tree.f90

    # Set iNodeMorton_I and iMortonNode_A indirect index arrays according to
    # 1. root node order
    # 2. Morton ordering for each root node

    #--------------------------------------------------------------------------
    print('h1')
    print(nNodeUsed)
    print('h2')

    #nNode = nRoot
    nNode_dict = {'nNode' : nRoot}
    #iNode = 0
    iNode_dict = {'iNode' : 0}
    #iMorton = 0
    iMorton_dict = {'iMorton' : 0}
    iNodeMorton_I = np.empty((nNodeUsed,))
    iNodeMorton_I[:] = Unset_
    iMortonNode_A = np.empty((nNodeUsed,))
    iMortonNode_A[:] = Unset_

    def order_children(iNode_tmp): # from recursive subroutine order_children in BATL_tree.f90

        # Recursively apply Morton ordering for nodes below a root block.
        # Store result into iNodeMorton_I and iMortonNode_A using the global
        # iMorton index.

        #--------------------------------------------------------------------------
        nNode_dict['nNode'] = np.maximum(nNode_dict['nNode'], iNode_tmp)

        if iTree_IA[FortEqu(Status_), FortEqu(iNode_tmp)] >= Used_:
            iMorton_dict['iMorton'] = iMorton_dict['iMorton'] + 1
            iNodeMorton_I[FortEqu(iMorton_dict['iMorton'])] = iNode_tmp
            #print('iNode_tmp = ' + str(iNode_tmp))
            iMortonNode_A[FortEqu(iNode_tmp)] = iMorton_dict['iMorton']
        else:
           for iChild in range(Child1_, ChildLast_+1):
                order_children(iTree_IA[FortEqu(iChild), FortEqu(iNode_tmp)])
        return None

    for kRoot in range(1, 1+nRoot_D[2]):
        for jRoot in range(1, 1+nRoot_D[1]):
            for iRoot in range(1, 1+nRoot_D[0]):
                #Root nodes are the first ones
                iNode_dict['iNode'] = iNode_dict['iNode'] + 1
                # All root nodes are handled as if they were first child
                order_children(iNode_dict['iNode'])


    #fortran: nNodeUsed = iMorton : yield UnboundLocalError
    assert(nNodeUsed == iMorton_dict['iMorton'])

    # Set min and max refinement levels
    nLevelMin = MaxLevel
    nLevelMax = 0
    for iMorton_tmp in range(1, nNodeUsed+1):
       iNode = iNodeMorton_I[FortEqu(iMorton_tmp)]
       iLevel = iTree_IA[FortEqu(Level_),FortEqu(iNode_dict['iNode'])]
       nLevelMin = np.minimum(iLevel, nLevelMin)
       nLevelMax = np.maximum(iLevel, nLevelMax)

    return nLevelMin, nLevelMax

# ############ from SWMF/GM/BATSRUS/srcBATL/BATL_size_orig.f90 
# # Maximum dimensionality of grid is 3 (cannot be modified)
# MaxDim = 3
# 
# # Indexes of AMR dimensions.
# # The magic formulas should be correct from 1 to nDimAmr.
# iDimAmrTmp_D = np.empty((MaxDim,))
# iDimAmrTmp_D = np.array([1 + (2-iRatio)*(3-jRatio), 6-iRatio-jRatio, 3 ], dtype=np.int8)
# if debug: print('iDimAmrTmp_D = ' + str(iDimAmrTmp_D))
# 
# iDimAmr_D = iDimAmrTmp_D[0:nDimAmr]
# ############
# if debug: print('iDimAmr_D = ' + str(iDimAmr_D))

# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975
def find_tree_node(CoordIn_D):
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

if debug: print(nNodeUsed)
#nLevelMin, nLevelMax = order_tree()
#print(nLevelMin, nLevelMax)
order_tree()
assert(False)

block_coords = []
position_mins = []
position_maxs = []
for iNode in blockused_:
    PositionMin_D, PositionMax_D, MaxIndex_D, BlockCoord_D = get_tree_position(iNode, returnall=True)

    block_coords.append(BlockCoord_D)
    position_mins.append(PositionMin_D)
    position_maxs.append(PositionMax_D)

block_coords = np.array(block_coords)
position_mins = np.array(position_mins)
position_maxs = np.array(position_maxs)

if debug:
    print(block_coords.shape)
    print(np.min(block_coords, axis=0))
    print(np.max(block_coords, axis=0))
    print(np.min(position_mins, axis=0))
    print(np.max(position_maxs, axis=0))


def get_physical_dimensions(iNode, returnCenterData=False):
    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]
    gridspacing = (256./8.)*0.5**iLevel

    PositionMin_D, PositionMax_D = get_tree_position(iNode)
    xmin = 256.*(PositionMin_D[0]) + -224.
    xmax = 256.*(PositionMax_D[0]) + -224.
    ymin = 256.*(PositionMin_D[1]) + -128.
    ymax = 256.*(PositionMax_D[1]) + -128.
    zmin = 256.*(PositionMin_D[2]) + -128.
    zmax = 256.*(PositionMax_D[2]) + -128.

    xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
    ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
    zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)

    xminmax = (xmin, xmax)
    yminmax = (ymin, ymax)
    zminmax = (zmin, zmax)

    if returnCenterData:
        return xlims, ylims, zlims, gridspacing
    else:
        return xminmax, yminmax, zminmax, gridspacing

counter = 0
for iNode in blockused_:
    xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(iNode, returnCenterData=False)

    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]
    if iLevel == 2:
        if debug:
            print('############\n')
            print(counter)
            #print((PositionMin_D, PositionMax_D))
            print(xminmax)
            print(yminmax)
            print(zminmax)
            print(gridspacing)
            if False:
                print('iNode = ' + str(iNode))

                print('Status_ --> ' + str(iTree_IA[FortEqu(Status_), FortEqu(iNode)]))
                print('Level_ --> ' + str(iTree_IA[FortEqu(Level_), FortEqu(iNode)]))
                print('Proc_ --> ' + str(iTree_IA[FortEqu(Proc_), FortEqu(iNode)]))
                print('Block_ --> ' + str(iTree_IA[FortEqu(Block_), FortEqu(iNode)]))
                print('MinLevel_ --> ' + str(iTree_IA[FortEqu(MinLevel_), FortEqu(iNode)]))
                print('MaxLevel_ --> ' + str(iTree_IA[FortEqu(MaxLevel_), FortEqu(iNode)]))
                print('Coord0_ --> ' + str(iTree_IA[FortEqu(Coord0_), FortEqu(iNode)]))
                print('Coord1_ --> ' + str(iTree_IA[FortEqu(Coord1_), FortEqu(iNode)]))
                print('Coord2_ --> ' + str(iTree_IA[FortEqu(Coord2_), FortEqu(iNode)]))
                print('Coord3_ --> ' + str(iTree_IA[FortEqu(Coord3_), FortEqu(iNode)]))
                print('CoordLast_ --> ' + str(iTree_IA[FortEqu(CoordLast_), FortEqu(iNode)]))
                print('Parent_ --> ' + str(iTree_IA[FortEqu(Parent_), FortEqu(iNode)]))
                print('Child0_ --> ' + str(iTree_IA[FortEqu(Child0_), FortEqu(iNode)]))
                print('Child1_ --> ' + str(iTree_IA[FortEqu(Child1_), FortEqu(iNode)]))
                print('ChildLast_ --> ' + str(iTree_IA[FortEqu(ChildLast_), FortEqu(iNode)]))

        counter += 1

    if counter>100: break

if __name__ == '__main__':
    import os
    import sys
    sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
    from config import conf
    from make_grid import make_axes, make_grid

    import spacepy.pybats.bats as bats
    import numpy as np

    mhd = bats.IdlFile(filetag + '.out')
    loaded_grid = np.column_stack([mhd['x'],mhd['y'],mhd['z']])
    npts = loaded_grid.shape[0]

    allgrids = []
    gridspacings = []
    counter = 0
    for iNode in blockused_:
        xlims, ylims, zlims, gridspacing = get_physical_dimensions(iNode, returnCenterData=True)

        axes = make_axes(xlims, ylims, zlims, gridspacing)
        grid = make_grid(axes)

        if counter<10:
            print(axes[0].shape,axes[1].shape,axes[2].shape)

        gridspacings.append(gridspacing)
        allgrids.append(grid)

        #ind = np.empty((nI*nJ*nK,))
        #for k in range(nI*nJ*nK):
        #    tr = loaded_grid == grid[k,:]
        #    ind[k] = np.where(tr[:,0] & tr[:,1] & tr[:,2])[0][0]

        counter += 1

    TreeGrid = np.array(allgrids)
    print('TreeGrid.shape == ' + str(TreeGrid.shape))
    reconstructed_grid = np.vstack(allgrids)
    print(np.all( reconstructed_grid.reshape(TreeGrid.shape) == TreeGrid ))
    print(reconstructed_grid)
    print(min(gridspacings))
    print(reconstructed_grid.shape)

    '''
    indecesReordered = np.empty((npts,))
    for k in range(npts):
        for ind in range(npts):
            bl = np.all(reconstructed_grid[k,:] == loaded_grid[ind,:])
            if bl:
                indecesReordered[k] = ind
                break

    print('checkpoint')
    assert(True)
    assert(np.all(  reconstructed_grid == loaded_grid[indecesReordered, :]  ))
    print('checkpoint2')
    '''

    #TreeVars = {}
    #TreeVars['p'] = np.empty((TreeGrid.shape[0],nI*nJ*nK))
    #for i in range(TreeGrid.shape[0]):
    #    TreeVars['p'][i,:] = mhd['p'][indexing[i]]
    #
    #print(TreeVars)

    reconstructed_grid = set([tuple(tup) for tup in list(reconstructed_grid)])
    loaded_grid = set([tuple(tup) for tup in list(loaded_grid)])

    print('\n\n\n\n######pleasebetrue####')
    print(reconstructed_grid == loaded_grid)
    print('\n\n')
    print(len(reconstructed_grid))


    #np.where(a[:,0] & a[:,1] & a[:,2])[0][0]

    '''from julia
    println(nDim)
    println(nInfo)
    println(nNode)
    println(iRatio_D)
    println(nRoot_D)
    println(size(iTree_IA))
    println(iTree_IA[:,1])
    println(iTree_IA[1,:])

    3
    18
    13161
    Int32[2, 2, 2]
    Int32[1, 1, 1]
    (18, 13161)
    Int32[-1, 0, -100, -100, 0, 30, 1, 1, 1, -100, 2, 3, 4, 5, 6, 7, 8, 9]
    '''

