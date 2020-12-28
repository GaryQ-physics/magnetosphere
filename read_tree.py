import numpy as np
import struct
##

# Number of cells per block in each direction.
# These values are set by the Config.pl script.
# Set 1 for ignored directions!
nI, nJ, nK = 4, 4, 4

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

# Named indexes of iTree_IA
Status_   =  1 - 1
Level_    =  2 - 1 # grid level
Proc_     =  3 - 1 # processor index
Block_    =  4 - 1 # block index
MinLevel_ =  5 - 1 # minimum level allowed
MaxLevel_ =  6 - 1 # maximum level allowed
Coord0_   =  6 - 1 # equal to Coord1_-1
Coord1_   =  7 - 1 # coordinate of node in 1st dimension
Coord2_   =  8 - 1 # coordinate of node in 2nd dimension
Coord3_   =  9 - 1 # coordinate of node in 3rd dimension
CoordLast_=  9 - 1 # Coord0_ + MaxDim (?)
Parent_   = 10 - 1 # Parent_ must be equal to Child0_
Child0_   = 10 - 1 #
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


filetag = "3d__var_3_e20031120-070000-000"



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
allblocks = [iTree_IA[Status_,iNode] for iNode in range(nNode)]
blockused_ = [iNode+1 for iNode in range(nNode) if iTree_IA[Status_,iNode] == Used_]
TEST = [iTree_IA[Proc_,iNode] for iNode in range(nNode) if iTree_IA[Status_,iNode] == Used_]

#print('#############\n\n')
#print(set(allblocks))
#print(set(TEST)==set(range(837))) for Proc_
#print(blockused_)





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

