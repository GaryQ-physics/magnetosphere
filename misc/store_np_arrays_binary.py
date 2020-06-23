"""
Demonstrates how to write a numpy array to a binary file and then how to read it.

It appears writing array to binary automatically flattens it. When you read in binary file that was written
    you get the flattened array 

for python2
"""

import numpy as np
import time

#large_array = np.linspace(0.,1e+5, 1e+10)

N=50000000
fname = 'np_array_test.bin'

to = time.time()
a = np.linspace(0., 1., N)
A = np.column_stack([a,a,a])
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

#print(A)
#print(A_load)
#print(A.flatten())

print( np.all(A.flatten() == A_load) )
