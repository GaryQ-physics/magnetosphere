import os
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from units_and_constants import phys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../swmf/')
import read_swmf_files as rswmf


def main(run, time, debug=True):
    if debug: print('filetag '+util.time2CDFfilename(run,time)[:-8])

    data, dTREE, ind, block2node, node2block = rswmf.read_all(util.time2CDFfilename(run,time)[:-8])
    nBlock, nI, nJ, nK = block2node.size, dTREE['nI'], dTREE['nJ'], dTREE['nK']
    block_data = {}
    for key in 'x y z jx jy jz bx by bz b1x b1y b1z'.split(' '):
        block_data[key] = data[key][ind].reshape((nBlock, nI, nJ, nK))

    gridspacing = np.zeros((nBlock, nI, nJ, nK), dtype=np.float32)
    epsilons = block_data['x'][:,1,0,0]-block_data['x'][:,0,0,0]
    gridspacing[:,:,:,:] = epsilons[:,None,None,None]
    block_data['gridspacing'] = gridspacing

    for vvar in ['j','b','b1']:
        if debug: print('computing partials for '+vvar)

        partials = np.zeros((nBlock, nI, nJ, nK, 3,3), dtype=np.float32)
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if i == 0 or i == nI-1:
                        partials[:,i,j,k, 0,0] = np.nan
                        partials[:,i,j,k, 0,1] = np.nan
                        partials[:,i,j,k, 0,2] = np.nan
                    else:
                        partials[:,i,j,k, 0,0] = (block_data[vvar+'x'][:, i+1, j  , k  ] - block_data[vvar+'x'][:, i-1, j  , k  ])/(2*epsilons)
                        partials[:,i,j,k, 0,1] = (block_data[vvar+'y'][:, i+1, j  , k  ] - block_data[vvar+'y'][:, i-1, j  , k  ])/(2*epsilons)
                        partials[:,i,j,k, 0,2] = (block_data[vvar+'z'][:, i+1, j  , k  ] - block_data[vvar+'z'][:, i-1, j  , k  ])/(2*epsilons)

                    if j == 0 or j == nJ-1:
                        partials[:,i,j,k, 1,0] = np.nan
                        partials[:,i,j,k, 1,1] = np.nan
                        partials[:,i,j,k, 1,2] = np.nan
                    else:
                        partials[:,i,j,k, 1,0] = (block_data[vvar+'x'][:, i  , j+1, k  ] - block_data[vvar+'x'][:, i  , j-1, k  ])/(2*epsilons)
                        partials[:,i,j,k, 1,1] = (block_data[vvar+'y'][:, i  , j+1, k  ] - block_data[vvar+'y'][:, i  , j-1, k  ])/(2*epsilons)
                        partials[:,i,j,k, 1,2] = (block_data[vvar+'z'][:, i  , j+1, k  ] - block_data[vvar+'z'][:, i  , j-1, k  ])/(2*epsilons)

                    if k == 0 or k == nK-1:
                        partials[:,i,j,k, 2,0] = np.nan
                        partials[:,i,j,k, 2,1] = np.nan
                        partials[:,i,j,k, 2,2] = np.nan
                    else:
                        partials[:,i,j,k, 2,0] = (block_data[vvar+'x'][:, i  , j  , k+1] - block_data[vvar+'x'][:, i  , j  , k-1])/(2*epsilons)
                        partials[:,i,j,k, 2,1] = (block_data[vvar+'y'][:, i  , j  , k+1] - block_data[vvar+'y'][:, i  , j  , k-1])/(2*epsilons)
                        partials[:,i,j,k, 2,2] = (block_data[vvar+'z'][:, i  , j  , k+1] - block_data[vvar+'z'][:, i  , j  , k-1])/(2*epsilons)

        if debug: print('computing vector derivatives for '+vvar)

        block_data['div_'+vvar] = partials[:,:,:,:,0,0]+partials[:,:,:,:,1,1]+partials[:,:,:,:,2,2]

        curl_tens = partials - partials.transpose((0,1,2,3,5,4))
        block_data['curl_'+vvar+'_x'] = curl_tens[:,:,:,:,1,2]
        block_data['curl_'+vvar+'_y'] = curl_tens[:,:,:,:,2,0]
        block_data['curl_'+vvar+'_z'] = curl_tens[:,:,:,:,0,1]

    unitmu0 = phys['mu0']*(phys['muA']/phys['m']**2)
    block_data['jRx'] = block_data['curl_b1_x']/unitmu0
    block_data['jRy'] = block_data['curl_b1_y']/unitmu0
    block_data['jRz'] = block_data['curl_b1_z']/unitmu0

    jRpartials = np.zeros((nBlock, nI, nJ, nK, 3,3), dtype=np.float32)
    for i in range(nI):
        for j in range(nJ):
            for k in range(nK):
                if i in [0, 1, nI-2, nI-1]:
                    jRpartials[:,i,j,k, 0,0] = np.nan
                    jRpartials[:,i,j,k, 0,1] = np.nan
                    jRpartials[:,i,j,k, 0,2] = np.nan
                else:
                    jRpartials[:,i,j,k, 0,0] = (block_data['jRx'][:, i+1, j  , k  ] - block_data['jRx'][:, i-1, j  , k  ])/(2*epsilons)
                    jRpartials[:,i,j,k, 0,1] = (block_data['jRy'][:, i+1, j  , k  ] - block_data['jRy'][:, i-1, j  , k  ])/(2*epsilons)
                    jRpartials[:,i,j,k, 0,2] = (block_data['jRz'][:, i+1, j  , k  ] - block_data['jRz'][:, i-1, j  , k  ])/(2*epsilons)

                if j in [0, 1, nJ-2, nJ-1]:
                    jRpartials[:,i,j,k, 1,0] = np.nan
                    jRpartials[:,i,j,k, 1,1] = np.nan
                    jRpartials[:,i,j,k, 1,2] = np.nan
                else:
                    jRpartials[:,i,j,k, 1,0] = (block_data['jRx'][:, i  , j+1, k  ] - block_data['jRx'][:, i  , j-1, k  ])/(2*epsilons)
                    jRpartials[:,i,j,k, 1,1] = (block_data['jRy'][:, i  , j+1, k  ] - block_data['jRy'][:, i  , j-1, k  ])/(2*epsilons)
                    jRpartials[:,i,j,k, 1,2] = (block_data['jRz'][:, i  , j+1, k  ] - block_data['jRz'][:, i  , j-1, k  ])/(2*epsilons)

                if k in [0, 1, nK-2, nK-1]:
                    jRpartials[:,i,j,k, 2,0] = np.nan
                    jRpartials[:,i,j,k, 2,1] = np.nan
                    jRpartials[:,i,j,k, 2,2] = np.nan
                else:
                    jRpartials[:,i,j,k, 2,0] = (block_data['jRx'][:, i  , j  , k+1] - block_data['jRx'][:, i  , j  , k-1])/(2*epsilons)
                    jRpartials[:,i,j,k, 2,1] = (block_data['jRy'][:, i  , j  , k+1] - block_data['jRy'][:, i  , j  , k-1])/(2*epsilons)
                    jRpartials[:,i,j,k, 2,2] = (block_data['jRz'][:, i  , j  , k+1] - block_data['jRz'][:, i  , j  , k-1])/(2*epsilons)

    block_data['div_jR'] = jRpartials[:,:,:,:,0,0]+jRpartials[:,:,:,:,1,1]+jRpartials[:,:,:,:,2,2]

    for key in block_data.keys():
        assert( block_data[key].shape == (nBlock,nI,nJ,nK) )

    for key in block_data.keys():
        block_data[key] = block_data[key].ravel()

    df_ends_as_nan = pd.DataFrame.from_dict(block_data)


    direct = conf[run+'_derived'] + 'derivatives/native_grid/'
    if not os.path.exists(direct): os.makedirs(direct)

    fname_df = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_df.pkl'%util.tpad(time, length=6)
    fname_meta = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_meta.txt'%util.tpad(time, length=6)

    meta = {'nBlock':nBlock, 'nI':nI, 'nJ':nJ, 'nK':nK}

    with open(fname_meta,'w') as handle:
        for key in meta.keys():
            handle.write('%s %d\n'%(key, meta[key]))
    if debug: print('saved ' + fname_meta)
    df_ends_as_nan.to_pickle(fname_df)
    if debug: print('saved ' + fname_df)

if __name__ == '__main__':
    ##################
    run = 'DIPTSUR2'
    debug = True

    if run == 'DIPTSUR2':
        time = (2019,9,2,6,30,0,0)
        #time = (2019,9,2,4,10,0,0)
    if run == 'IMP10_RUN_SAMPLE':
        time = (2019,9,2,7,0,0,0)
    if run == 'TESTANALYTIC':
        time = (2000,1,1,0,10,0,0)
    if run == 'LUHMANN1979':
        time = (2000,1,1,0,0,0,0)
    ##################

    main(run, time, debug=debug)
