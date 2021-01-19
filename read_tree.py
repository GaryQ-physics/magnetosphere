import numpy as np
import struct

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

# from SWMF/GM/BATSRUS/srcBATL_/BATL_tree.f90
MaxCoord_I =  [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 
                4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
                268435456, 536870912, 1073741824 ]
assert(len(MaxCoord_I) == MaxLevel+1)

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


filetag = "/home/gary/Batsrus.jl-master/3d__var_3_e20031120-070000-000"



## Loading AMR tree

# directly read bytes
f = open(filetag+".tree", 'rb')
filebytes = np.fromfile(f, dtype=np.int32)
f.close()
print(filebytes[:20])

if False:
	# use scipy FortranFile
	from scipy.io import FortranFile
	ff = FortranFile(filetag+".tree", 'r')
	if False:
		nDim, nInfo, nNode = ff.read_reals(dtype=np.int32)
		iRatio_D = ff.read_reals(dtype=np.int32) # Array of refinement ratios
		nRoot_D = ff.read_reals(dtype=np.int32) # The number of root nodes in all dimension
		iTree_IA = ff.read_reals(dtype=np.int32).reshape((nInfo,nNode), order='fortran')
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


'''WRONG
headerInt = filebytes[0]
nDim = filebytes[1]
nInfo = filebytes[2]
nNode = filebytes[3]
iRatio_D = filebytes[4:4+nDim] # Array of refinement ratios
nRoot_D = filebytes[4+nDim:4+2*nDim]  # The number of root nodes in all dimension
iTree_IA = filebytes[4+2*nDim:4+2*nDim+nInfo*nNode].reshape((nInfo,nNode))
'''

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

MinLevels = [iTree_IA[FortEqu(MinLevel_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
MaxLevels = [iTree_IA[FortEqu(MaxLevel_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Blocks = [iTree_IA[FortEqu(Block_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Levels = [iTree_IA[FortEqu(Level_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]
Procs = [iTree_IA[FortEqu(Proc_), FortEqu(iNode)] for iNode in range(1, nNode+1) if iTree_IA[FortEqu(Status_), FortEqu(iNode)] == Used_]

print((256./8.)*0.5**max(Levels))
print((256./8.)*0.5**min(Levels))


print('#############\n\n')
print(set(MinLevels))
print(set(MaxLevels))
print('#############\n\n')

#print('#############\n\n')
#print(set(allblocks))
#print(set(TEST)==set(range(837))) for Proc_
#print(blockused_)

# from SWMF/GM/BATSRUS/srcBATL_/BATL_tree.f90 line 951 (they dont return MaxIndex_D
def get_tree_position(iNode, returnall=False):
    # Calculate normalized position of the edges of node inode.
    # Zero is at the minimum boundary of the grid, one is at the max boundary

    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]

    # For non-AMR directions MaxIndex_D = nRoot_D
    # For AMR     directions MaxIndex_D = nRoot_D*MaxCoord_I(iLevel)
    MaxIndex_D = ((MaxCoord_I[FortEqu(iLevel)]-1)*(iRatio_D-1) + 1)*nRoot_D

    # Convert to real by adding -1.0 or 0.0 for the two edges, respectively
    PositionMin_D = (iTree_IA[FortEqu(Coord1_):FortEqu(CoordLast_)+1,FortEqu(iNode)] - 1.0)/MaxIndex_D
    PositionMax_D = (iTree_IA[FortEqu(Coord1_):FortEqu(CoordLast_)+1,FortEqu(iNode)] + 0.0)/MaxIndex_D

    if returnall:
        return PositionMin_D, PositionMax_D, MaxIndex_D
    else:
        return PositionMin_D, PositionMax_D


counter = 0
for iNode in blockused_:
    PositionMin_D, PositionMax_D = get_tree_position(iNode)
    xmin = 256.*PositionMin_D[0] + -224.
    xmax = 256.*PositionMax_D[0] + -224.
    ymin = 256.*PositionMin_D[1] + -128.
    ymax = 256.*PositionMax_D[1] + -128.
    zmin = 256.*PositionMin_D[2] + -128.
    zmax = 256.*PositionMax_D[2] + -128.
    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]
    gridspacing = (256./8.)*0.5**iLevel
    xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
    ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
    zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)

    if iLevel == 2:
        print('############\n')
        print(counter)
        print((PositionMin_D, PositionMax_D))
        print(xlims)
        print(ylims)
        print(zlims)
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


assert(False)


from config import conf
from make_grid import make_axes, make_grid
counter = 0
allgrid = []
gridspacings = []
for iNode in blockused_:
    PositionMin_D, PositionMax_D = get_tree_position(iNode)
    xmin = 256.*PositionMin_D[0] + -224.
    xmax = 256.*PositionMax_D[0] + -224.
    ymin = 256.*PositionMin_D[1] + -128.
    ymax = 256.*PositionMax_D[1] + -128.
    zmin = 256.*PositionMin_D[1] + -128.
    zmax = 256.*PositionMax_D[1] + -128.
    iLevel = iTree_IA[FortEqu(Level_), FortEqu(iNode)]
    gridspacing = (256./8.)*0.5**iLevel
    xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
    ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
    zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)

    axes = make_axes(xlims, ylims, zlims, gridspacing)
    grid = make_grid(axes)

    if counter<100:
        print(axes[0].shape,axes[1].shape,axes[2].shape)

    gridspacings.append(gridspacing)
    allgrid.append(grid)

    counter += 1

allgrid = np.vstack(allgrid)
print(allgrid)
print(min(gridspacings))
print(allgrid.shape)


#print(iTree_IA[FortEqu(Level_), FortEqu(iNode)])
#print(iTree_IA[FortEqu(Coord1_):FortEqu(CoordLast_)+1,FortEqu(iNode)])


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

