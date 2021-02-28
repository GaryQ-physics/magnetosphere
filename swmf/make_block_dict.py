import os
import sys
sys.path.append('/home/gary/magnetosphere/physics/')
from make_grid import make_axes, make_grid
import numpy as np
import spacepy.pybats.bats as bats

fname = "../../Batsrus.jl-master/3d__var_1_n00002500"

data = bats.Bats2d(fname + ".out")
points = np.column_stack([data['x'],data['y'],data['z']])
ind = np.load(fname+'.out-indices_for_reordering.npy')
ind
TreeGrid = points[ind,:].reshape((11516, 512, 3))

TreeGrid.shape
np.all(TreeGrid[:,:, :].flatten() == points[ind,:].flatten())
np.all(TreeGrid[:,:, 0].flatten() == points[ind,0].flatten())
np.all(TreeGrid[:,:, 0].flatten() == points[ind,0])
np.all(TreeGrid[:,:,0].flatten() == points[ind,0])
np.all(TreeGrid[:,:,0].flatten() == data['x'][ind])
np.all(TreeGrid[:,:,0] == data['x'][ind].reshape((11516,512)))
np.all(TreeGrid[:,:,1] == data['y'][ind].reshape((11516,512)))
block_data={}

block_data['x'] = data['x'][ind].reshape((11516,8,8,8))
block_data['y'] = data['y'][ind].reshape((11516,8,8,8))
block_data['z'] = data['z'][ind].reshape((11516,8,8,8))
block_data['bx'] = data['bx'][ind].reshape((11516,8,8,8))
block_data['by'] = data['by'][ind].reshape((11516,8,8,8))
block_data['bz'] = data['bz'][ind].reshape((11516,8,8,8))
block_data['p'] = data['p'][ind].reshape((11516,8,8,8))

block_data['b1x'] = data['b1x'][ind].reshape((11516,8,8,8))
block_data['b1y'] = data['b1y'][ind].reshape((11516,8,8,8))
block_data['b1z'] = data['b1z'][ind].reshape((11516,8,8,8))

block_data['jx'] = data['jx'][ind].reshape((11516,8,8,8))
block_data['jy'] = data['jy'][ind].reshape((11516,8,8,8))
block_data['jz'] = data['jz'][ind].reshape((11516,8,8,8))

block_data['rho'] = data['rho'][ind].reshape((11516,8,8,8))
