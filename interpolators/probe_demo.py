import numpy as np
from probe import probe, GetRunData

'''
fname_templ = '/home/gary/temp/'+'3d__var_3_e20031120-070000-000'
d=probe(fname_templ+'.vtu', (-10, 0, 0), var='p', library='vtk')
print(type(d))
print(d)

d = GetRunData('SCARR5', (2003, 11, 20, 7, 0, 0), (-10, 0, 0), 'p')
print(type(d))
print(d)

d = GetRunData('SCARR5', (2003, 11, 20, 7, 0, 0), (-10, 0, 0), 'b')
print(type(d))
print(d)
'''
X = np.array([[-5, 0, 0], [-10, 0, 0]])
#j = GetRunData('SCARR5', (2003, 11, 20, 7, 0, 0), (-10, 0, 0), 'j')

print(X.shape)

print('\n\n\n')

d_kam = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), (-10, 0, 0), 'p', library='kameleon')
d_vtk = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), (-10, 0, 0), 'p', library='vtk')
print(type(d_kam))
print(type(d_vtk))
print(d_kam)
print(d_vtk)

d_kam = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), (-10, 0, 0), 'b1', library='kameleon')
d_vtk = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), (-10, 0, 0), 'b1', library='vtk')
print(type(d_kam))
print(type(d_vtk))
print(d_kam)
print(d_vtk)

d_kam = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), X, 'p', library='kameleon')
d_vtk = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), X, 'p', library='vtk')
print(type(d_kam))
print(type(d_vtk))
print(d_kam)
print(d_vtk)

d_kam = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), X, 'b1', library='kameleon')
d_vtk = GetRunData('DIPTSUR2', (2019,9,2,6,30,0,0), X, 'b1', library='vtk')
print(type(d_kam))
print(type(d_vtk))
print(d_kam)
print(d_vtk)
