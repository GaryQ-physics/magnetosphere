import numpy as np
from probe import probe, probe_vect

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var='p')
print(d)

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var=['p'])
print(d)

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var=['p', 'bx'])
print(d)

'''
jx = probe((2003, 11, 20, 7, 0, 0), (-5,0,0), var='jx', debug=True)
jy = probe((2003, 11, 20, 7, 0, 0), (-5,0,0), var='jy', debug=True)
jz = probe((2003, 11, 20, 7, 0, 0), (-5,0,0), var='jz', debug=True)
j = probe((2003, 11, 20, 7, 0, 0), (-5,0,0), var='j', debug=True)
print(jx, jy, jz)
print(j)
print(np.sqrt(jx**2 + jy**2 + jz**2))
'''

X = np.array([[-5, 0, 0], [-10, 0, 0]])

jy = probe((2003, 11, 20, 7, 0, 0), X, var='jy')

j_dict = probe((2003, 11, 20, 7, 0, 0), X, var=['jx', 'jy', 'jz'])

#j = np.column_stack([j_dict['jx'], j_dict['jy'], j_dict['jz']])
j = probe_vect((2003, 11, 20, 7, 0, 0), X, 'j')

print(j_dict)
print(X.shape)
print(j.shape)
print(j)
