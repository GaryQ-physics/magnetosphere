"""
Demonstrates how to write a numpy array to a binary file and then how to read it.

Writing array to binary automatically flattens it. When you read in binary file
that was written you get the flattened array back.

for python2

Typical out:
    N = 50000000
    time to generate Nx3 array = 0.85 sec
    time to write Nx3 array to binary file = T 
    time to read binary file back to numpy array = 0.26 sec
    ('shape of resulting np array is', (150000000,))
    ('np.all(A.flatten() == A_load) -->', True)
NOTE: where here T depends strongly on wether the file np_array_test.bin already exists.
    if it doesn't already exist, then it takes about 1 to 3 seconds
    if it does and thus needs to be rewritten, it takes about 15 to 18 sec

"""
import os
import time
import tempfile
import numpy as np

N = 50000000
fname = tempfile.gettempdir() + '/np_array_test.bin'

print('N = {0:d}'.format(N))

to = time.time()
a = np.linspace(0., 1., N)
A = np.column_stack([a, a, a])
tf = time.time()
print('time to generate Nx3 array = {0:.5f} sec'.format(tf-to))

#print(A.shape)

to = time.time()

if os.path.exists(fname):
    os.remove(fname)

if False:
    f = open(fname, 'w')
    f.write(A.tobytes())
    f.close()
else:
    A.tofile(fname)

tf = time.time()
print('time to write Nx3 array to binary file = {0:.5f} sec'.format(tf-to))

to = time.time()
A_load = np.fromfile(fname)
tf = time.time()
print('time to read binary file back to numpy array = {0:.5f} sec'.format(tf-to))


print('shape of resulting np array is', A_load.shape)
print('np.all(A.flatten() == A_load) -->', np.all(A.flatten() == A_load) )
