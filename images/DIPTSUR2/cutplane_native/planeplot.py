import numpy as np
import matplotlib.pyplot as plt

planedata = np.load('20190902T041000_x_z_divB1_normB1.npy')
print(planedata.shape)

x = planedata[0,:]
z = planedata[1,:]
divB1 = planedata[2,:]
normB1 = planedata[3,:]

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
plt.title('y = 3/32 R_E plane')
plt.xlabel('x [R_E]')
plt.ylabel('z [R_E]')
plt.savefig('div_b1_cutplane.png')
#plt.show()
plt.clf()

plt.figure(figsize=(8,6))
plt.scatter(x,z,c=normB1)
plt.colorbar(label='norm_b1 [$nT$]')
plt.title('y = 3/32 R_E plane')
plt.xlabel('x [R_E]')
plt.ylabel('z [R_E]')
plt.savefig('norm_b1_cutplane.png')
#plt.show()
plt.clf()
