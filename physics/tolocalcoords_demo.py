import biot_savart_kameleon_interpolated_grid as bsk
import numpy as np

time = (2019,9,2,6,30,0)

dB = np.array([[1.,2,3],[4,5,6]])
dB_loc = bsk.toMAGLocalComponents(time, 11., 147., dB)

U1 = np.array([ 0.14593682, 0.98929282, 0.0014703])
U2 = np.array([-0.19864568,  0.02784747,  0.97967567])
U3 = np.array([0.96914516, -0.14326282,  0.20058272])
M = np.column_stack([U2,U1,-U3])
print(M)

Minv = np.linalg.inv(M)

print(np.matmul(Minv,dB[0,:]))
print(np.einsum('ij,j', Minv,dB[0,:]))

print(np.matmul(Minv,dB[1,:]))
print(np.einsum('ij,j', Minv,dB[1,:]))

print(dB_loc)
