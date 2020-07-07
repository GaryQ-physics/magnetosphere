import numpy as np
from probe import probe

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var='p')
print(d)

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var=['p'])
print(d)

d = probe((2003, 11, 20, 7, 0, 0), (-10, 0, 0), var=['p', 'bx'], dictionary=True)
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

jy1 = probe((2003, 11, 20, 7, 0, 0), [-5, 0, 0], var='jy')
jy2 = probe((2003, 11, 20, 7, 0, 0), [-10, 0, 0], var='jy')

print(jy1)
print(jy2)

X = np.array([[-5, 0, 0], [-10, 0, 0]])
jy = probe((2003, 11, 20, 7, 0, 0), X, var='jy')

print(jy.shape)
print(jy)

j = probe((2003, 11, 20, 7, 0, 0), X, var=['jx', 'jy', 'jz'])

print(X.shape)
print(j.shape)
print(j)


#j = np.column_stack([j_dict['jx'], j_dict['jy'], j_dict['jz']])
