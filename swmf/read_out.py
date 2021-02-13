'''tested with python2.7 and python3.7
with python2.7, works for either usesp True or False
with python3.7, only works for either usesp True 
    (there is no fortranfile library installed, however
     pip install fortranfile was tried, the install worked, but the
     script crashed at runtime)

with usesp = True, the header is stored as an int8 array, rather than a string
'''
import numpy as np

usesp = True
nDim = 3
filetag = "3d__var_3_e20031120-070000-000"
filetag = "../../Batsrus.jl-master/3d__var_1_n00002500"

# directly read bytes
f = open(filetag+".out", 'rb')
filebytes = np.fromfile(f, dtype=np.int32)
f.close()
print(filebytes[:20])

if usesp:
    # use scipy FortranFile
    from scipy.io import FortranFile
    ff = FortranFile(filetag+".out", 'r')
    head = ff.read_ints(dtype=np.int8)
    nStep, Time, nDimOut, nParam, nVar = ff.read_ints(dtype=np.int32)
    n_D = ff.read_ints(dtype=np.int32)
else:
    # use fortranfile
    from fortranfile import FortranFile
    ff = FortranFile(filetag+".out")
    head = ff.readString()
    nStep, Time, nDimOut, nParam, nVar = ff.readInts()  # Time should be Real4_ ???
    n_D = ff.readInts()
    
     



print(head)
print((nStep, Time, nDimOut, nParam, nVar))
print(n_D)

if usesp:
    A1 = ff.read_ints(dtype=np.int32)
    A2 = ff.read_ints(dtype=np.int32)
else:
    A1 = ff.readInts()
    A2 = ff.readInts()

print(A1.shape)
print(A1)
print(A2.shape)
print(A2)

A = None
for i in range(50):
    if usesp:
        A = ff.read_reals(dtype=np.float32)
    else:
        A = ff.readReals()
    print(A.shape)
