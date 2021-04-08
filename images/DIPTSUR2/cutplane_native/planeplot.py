import numpy as np
import matplotlib.pyplot as plt

planedata = np.load('y=3_16_20190902T063000_x_z_divB1_normB1_normcurlB1.npy')
title = 'y = 3/16 R_E plane  at 6:30'
name = 'y=3_16_at_6_30'

print(planedata.shape)

x = planedata[0,:]
z = planedata[1,:]
divB1 = planedata[2,:]
normB1 = planedata[3,:]
normcurlB1 = planedata[4,:]

rcut=2.
tr = 0.09375**2 + x**2 + z**2 <= rcut**2

x[tr] = np.nan
z[tr] = np.nan
divB1[tr] = np.nan
divB1[tr] = np.nan

print(np.count_nonzero(np.isnan(x)))
print(np.count_nonzero(np.isnan(z)))
print(np.count_nonzero(np.isnan(divB1)))
print(np.count_nonzero(np.isnan(normB1)))

#x = np.array([1.,1.,2.,2.])
#z = np.array([1.,2.,1.,2.])
#divB1 = np.array([1.,2.,3.,np.nan])

plt.figure(figsize=(8,6))
plt.scatter(x,z,c=divB1)
plt.colorbar(label='div_b1 [$nT R_E^{-1}$]')
plt.title(title)
plt.xlabel('x [R_E]')
plt.ylabel('z [R_E]')
plt.savefig(name+'-div_b1.png')
#plt.show()
#plt.clf()

plt.figure(figsize=(8,6))
plt.scatter(x,z,c=normB1)
plt.colorbar(label='norm_b1 [$nT$]')
plt.title(title)
plt.xlabel('x [R_E]')
plt.ylabel('z [R_E]')
plt.savefig(name+'-norm_b1.png')
#plt.show()
#plt.clf()

plt.figure(figsize=(8,6))
plt.scatter(x,z,c=normcurlB1)
plt.colorbar(label='norm_curl_b1 [$nT R_E^{-1}$]')
plt.title(title)
plt.xlabel('x [R_E]')
plt.ylabel('z [R_E]')
plt.savefig(name+'-norm_curl_b1.png')
#plt.show()
#plt.clf()
