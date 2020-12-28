import numpy as np

nDim = 3
filetag = "3d__var_3_e20031120-070000-000"

# directly read bytes
f = open(filetag+".out", 'rb')
filebytes = np.fromfile(f, dtype=np.int32)
f.close()
print(filebytes[:20])

if False:
	# use scipy FortranFile
	from scipy.io import FortranFile
	ff = FortranFile(filetag+".out", 'r')
	head = ff.read_record(dtype=str)
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
