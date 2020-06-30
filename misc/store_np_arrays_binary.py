"""
Demonstrates how to write a numpy array to a binary file and then how to read it.

Writing array to binary automatically flattens it. When you read in binary file
that was written you get the flattened array back.

for python2
"""

import time
import numpy as np

N = 50000
fname = '/tmp/np_array_test.bin'

to = time.time()
a = np.linspace(0., 1., N)
A = np.column_stack([a, a, a])
tf = time.time()
print(tf-to) 

print(A.shape)

to = time.time()
f = open(fname, 'w')
f.write(A.tobytes())
f.close()
tf = time.time()
print(tf-to) 

to = time.time()
A_load = np.fromfile(fname)
tf = time.time()
print(tf-to) 

print(A_load.shape)
print( np.all(A.flatten() == A_load) )
