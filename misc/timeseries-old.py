'''man
 -s : summaries
 -c : coulomb
 -b : biot savart
 -p : compute in parallel
 -v : verbose, print details ect
 --nfiles : set equal to n_files

eg:
python timeseries.py -scb --nfiles=4
'''
import os
import sys

n_times =  None
do_summary = False
do_coulomb = False
do_COULOMB = False
do_biotsavart = False
do_BIOTSAVART = False
do_inmag = False
para = False
verbose = False
for arg in sys.argv:
    if arg[0:2] == '--':
        if '--nfiles' == arg[0:8]:
            n_times = int(arg[9:])
        else:
            raise ValueError ('unknown argument '+arg)
    elif arg[0] == '-':
        for char in arg[1:]:
            if char=='s':
                do_summary = True
            if char=='c':
                do_coulomb = True
            if char=='C':
                do_COULOMB = True
            if char=='b':
                do_biotsavart = True
            if char=='B':
                do_BIOTSAVART = True
            if char=='m':
                do_inmag= True
            if char=='p':
                para = True
            if char=='v':
                verbose = True

import numpy as np
import pandas as pd
import datetime
from config import conf
import util
import read_swmf_files as rswmf
from magnetometers import GetMagnetometerLocation
from timeseries_summarize2 import get_summary
from named_var_indexes import index2str, nVarTot,_b1x,_b1y,_b1z
from coulomb import B_biotsavart, B_coulomb

run = 'DIPTSUR2'

rcut = util.get_rCurrents(run)
direct = conf[run+'_derived'] + 'timeseries/'
if not os.path.exists(direct): os.makedirs(direct)

# wrapper function, to run all functions desired by command line argument,
#    while only loading the .out once per time slice. 
def wrap(time):
    NeededArray = rswmf.get_needed_array(run,time)

    if do_summary:
        summary_arr = get_summary(NeededArray, rcut)

        if not os.path.exists(direct+'do_summary/'): os.makedirs(direct+'do_summary/')
        outname = direct+'do_summary/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_summary_arr.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, summary_arr)

    if do_coulomb:
        Bcl = B_coulomb( np.zeros(3), NeededArray, rcut=rcut)

        if not os.path.exists(direct+'do_coulomb/'): os.makedirs(direct+'do_coulomb/')
        outname = direct+'do_coulomb/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_Bcl.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, Bcl)

    if do_biotsavart:
        Bbs = B_biotsavart( np.zeros(3), NeededArray, rcut=rcut)

        if not os.path.exists(direct+'do_biotsavart/'): os.makedirs(direct+'do_biotsavart/')
        outname = direct+'do_biotsavart/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_Bbs.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, Bbs)

    if do_COULOMB:
        colaba = np.array(GetMagnetometerLocation('colaba', time, 'GSM', 'car'))
        BCL = B_coulomb( colaba, NeededArray, rcut=rcut)

        if not os.path.exists(direct+'do_COULOMB/'): os.makedirs(direct+'do_COULOMB/')
        outname = direct+'do_COULOMB/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_BCL.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, BCL)

    if do_BIOTSAVART:
        colaba = np.array(GetMagnetometerLocation('colaba', time, 'GSM', 'car'))
        BBS = B_biotsavart( colaba, NeededArray, rcut=rcut)

        if not os.path.exists(direct+'do_BIOTSAVART/'): os.makedirs(direct+'do_BIOTSAVART/')
        outname = direct+'do_BIOTSAVART/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_BBS.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, BBS)

    if do_inmag:
        point = np.array([0.71875,0.09375,-3.71875])
        assert(point[0] == NeededArray[0,9477, 3,1,4])
        assert(point[1] == NeededArray[1,9477, 3,1,4])
        assert(point[2] == NeededArray[2,9477, 3,1,4])

        b1x_nat = NeededArray[_b1x,9477, 3,1,4]
        b1y_nat = NeededArray[_b1y,9477, 3,1,4]
        b1z_nat = NeededArray[_b1z,9477, 3,1,4]

        Bcl = B_coulomb( point, NeededArray, rcut=rcut)
        Bbs = B_biotsavart( point, NeededArray, rcut=rcut)
        integrals = np.hstack([Bcl, Bbs, b1x_nat,b1y_nat,b1z_nat])

        if not os.path.exists(direct+'do_inmag/'): os.makedirs(direct+'do_inmag/')
        outname = direct+'do_inmag/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_integrals.npy'%util.tpad(time, length=6)
        if verbose: print('saving '+outname)
        np.save(outname, integrals)   


# loop through each time slice, in parallel or in serial, and execute wrapper
times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:1*n_times:1]

if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is None:
        assert(False)
    num_cores = min(num_cores, len(times), 20)
    print('Parallel processing {0:d} time slices using {1:d} cores'\
          .format(len(times), num_cores))
    Parallel(n_jobs=num_cores)(\
              delayed(wrap)(time) for time in list(times))
else:
    for time in list(times):
        wrap(time)

# stitch the written files together
if do_summary:
    epsilons = [ 0.0625,
                 0.1250,
                 0.2500,
                 0.5000,
                 1.0000,
                 2.0000,
                 4.0000,
                 8.0000 ]
    _count = 0
    _mean  = 1
    _sndMo = 2
    _std   = 3
    _min   = 4
    _max   = 5 #!!! warning: duplicated

    list_summaries = []
    dtimes = []
    for time in times:
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_summary/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_summary_arr.npy'%util.tpad(time, length=6)

        list_summaries.append(np.load(outname))

    multiind = pd.MultiIndex.from_product([epsilons, ['count', 'mean', 'sndMo', 'std', 'min', 'max']])
    for i_var in range(nVarTot):
        summaryDF = pd.DataFrame(index=dtimes, columns=multiind)

        for i_eps in range(len(epsilons)):
            for i_t in range(len(times)):
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'count')] = list_summaries[i_t][i_var,i_eps, _count]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'mean' )] = list_summaries[i_t][i_var,i_eps, _mean ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'sndMo')] = list_summaries[i_t][i_var,i_eps, _sndMo]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'std'  )] = list_summaries[i_t][i_var,i_eps, _std  ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'min'  )] = list_summaries[i_t][i_var,i_eps, _min  ]
                summaryDF.loc[dtimes[i_t], (epsilons[i_eps],'max'  )] = list_summaries[i_t][i_var,i_eps, _max  ]

        summaryDF.to_pickle(direct+'summaryDF_%s.pkl'%(index2str[i_var]))

if do_coulomb:
    dtimes = []
    Bcls = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_coulomb/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_Bcl.npy'%util.tpad(time, length=6)
        Bcls.append(np.load(outname))

    df = pd.DataFrame(data=Bcls, columns = ['B_coulomb_x','B_coulomb_y','B_coulomb_z'],
                        index=dtimes)
    df.to_pickle(direct+'df_coulomb_integral_for_earth_center.pkl')

if do_biotsavart:
    dtimes = []
    Bbss = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_biotsavart/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_Bbs.npy'%util.tpad(time, length=6)
        Bbss.append(np.load(outname))

    df = pd.DataFrame(data=Bbss, columns = ['B_biotsavart_x','B_biotsavart_y','B_biotsavart_z'],
                        index=dtimes)
    df.to_pickle(direct+'df_biotsavart_integral_for_earth_center.pkl')

if do_COULOMB:
    dtimes = []
    BCLs = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_COULOMB/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_BCL.npy'%util.tpad(time, length=6)
        BCLs.append(np.load(outname))

    df = pd.DataFrame(data=BCLs, columns = ['B_coulomb_x','B_coulomb_y','B_coulomb_z'],
                        index=dtimes)
    df.to_pickle(direct+'df_coulomb_integral_for_colaba.pkl')

if do_BIOTSAVART:
    dtimes = []
    BBSs = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_BIOTSAVART/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_BBS.npy'%util.tpad(time, length=6)
        BBSs.append(np.load(outname))

    df = pd.DataFrame(data=BBSs, columns = ['B_biotsavart_x','B_biotsavart_y','B_biotsavart_z'],
                        index=dtimes)
    df.to_pickle(direct+'df_biotsavart_integral_for_colaba.pkl')

if do_inmag:
    dtimes = []
    integrals = []
    for time in list(times):
        dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
        outname = direct+'do_inmag/' \
                        + '%.2d%.2d%.2dT%.2d%.2d%.2d_integrals.npy'%util.tpad(time, length=6)
        integrals.append(np.load(outname))

    df = pd.DataFrame(data=integrals, columns = ['B_coulomb_x','B_coulomb_y','B_coulomb_z',
                                            'B_biotsavart_x','B_biotsavart_y','B_biotsavart_z',
                                            'b1x','b1y','b1z'],
                        index=dtimes)
    df.to_pickle(direct+'df_integrals_in_magnetosphere.pkl')
