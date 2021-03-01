import numpy as np
import read_swmf_files as rswmf

filetag = "/home/gary/temp/3d__var_3_e20031120-070000-000"

dTREE = rswmf.read_tree_file(filetag)
rswmf.test(dTREE) # should have no errors

points = 10.*np.random.rand(33).reshape((11,3))
for i in range(11):
    iNode = rswmf.find_tree_node(points[i,:], dTREE)
    dims = rswmf.get_physical_dimensions(iNode, dTREE)
    print('\npoint = %s\nat iNode = %d\nwith dims = %s\n'%(str(points[i,:]),iNode,str(dims)))



import spacepy.pybats.bats as bats
import read_swmf_files as rswmf
data = bats.Bats2d(filetag + ".out")
points = np.column_stack([data['x'],data['y'],data['z']])
print('ready\n')

for i in [230000,420000,620000]:
    point = points[i,:]
    actual = data['p'][i]
    interp = rswmf.interpolate(filetag, point,var='p')
    print('\npoint = %s\nwith value = %s\ninterpolated to = %s\nareEqual: %s\n'%(
            str(point),
            str(actual),
            str(interp),
            str(actual==interp)
            ))




def test(dTREE):
    assert(dTREE['iTree_IA'][FortEqu(Level_), FortEqu(iRootNode)] == 0)
    xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(iRootNode, dTREE, returnCenters=False)
    assert(xminmax == (dTREE['xGlobalMin'], dTREE['xGlobalMax']))
    assert(yminmax == (dTREE['yGlobalMin'], dTREE['yGlobalMax']))
    assert(zminmax == (dTREE['zGlobalMin'], dTREE['zGlobalMax']))

