import os
import sys
import tempfile
import resource

if os.path.exists('/Users/robertweigel/'):
    base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'

if os.path.exists('/Users/weigel/'):
    base = '/Users/weigel/git/magnetosphere/'

elif os.path.exists('/home/weigel/') and False:
    base = '/home/weigel/git/magnetosphere/'

elif os.path.exists('/home/gary/'):
    base = '/home/gary/magnetosphere/'
    storage = base +'data/'
    #https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
    soft, hard = int(13*2**30), int(13*2**30)
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))

elif os.path.exists('/home/gquaresi/'):
    base = '/home/gquaresi/magnetosphere/'
    storage = '/media/sunspot/git-data/sblake/'
    #https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
    soft, hard = 90*2**30, 90*2**30
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
else:
    assert(False)


defined_magnetometers = {
    'YKC'      : ('GEO','sph', 1.,  62.480 ,  245.518 ), # radius, lat, lon   in GEO
    'FRN'      : ('GEO','sph', 1.,  37.0913, -119.7193),
    'FUR'      : ('GEO','sph', 1.,  48.17  ,  11.28   ),
    'colaba'   : ('GEO','sph', 1.,  18.907 ,  72.815  ),
    'GMpoint1' : ('GSM','car', 1.,  18.907 ,  72.815  ),

    }


conf = {
        'storage': storage,
        'base': base,
        'SWPC_raw': storage + 'SWPC_SWMF_052811_2/raw_output/',
        'SWPC_cdf_path': storage + 'SWPC_SWMF_052811_2/GM_CDF/',

        'SWPC_cdf': storage + 'SWPC_SWMF_052811_2/GM_CDF/',
        'SWPC_derived': storage + 'SWPC_SWMF_052811_2-derived/',

        'SCARR5_cdf': storage + 'SCARR5/GM/IO2/',
        'SCARR5_derived': storage + 'SCARR5-derived/',
        'SCARR5_magfile': storage + 'SCARR5/MAG_FILES/',

        'SCARR1_cdf': storage + 'SCARR1/tmp/',
        'SCARR1_magfile': storage + 'SCARR1/MAG_FILES/',
        'SCARR1_derived': storage + 'SCARR1-derived/',

        'CARR_IMPULSE_cdf': storage + 'CARR_IMPULSE_3D/GM/IO2/',
        'CARR_IMPULSE_derived': storage + 'CARR_IMPULSE_3D-derived/',

        'DIPTSUR2_cdf' : storage + 'DIPTSUR2/GM/IO2/',
        'DIPTSUR2_derived' : storage + 'DIPTSUR2-derived/',

        'DIPTSUR2_4HOUR_cdf' : storage + 'DIPTSUR2_4HOUR/GM/IO2/',
        'DIPTSUR2_4HOUR_derived' : storage + 'DIPTSUR2_4HOUR-derived/',

        'UNALT_DIPOLE_cdf' : storage + 'UNALT_DIPOLE/GM/IO2/',
        'UNALT_DIPOLE_derived' : storage + 'UNALT_DIPOLE-derived/',

        'UNALT_DIPOLE4_cdf' : storage + 'UNALT_DIPOLE4/GM/IO2/',
        'UNALT_DIPOLE4_derived' : storage + 'UNALT_DIPOLE4-derived/',

        'LUHMANN1979_cdf' : storage + 'LUHMANN1979/',
        'LUHMANN1979_derived' : storage + 'LUHMANN1979-derived/',

        'mag_server_url': 'http://mag.gmu.edu/git-data/sblake/',

        'SWPC_raw': storage + 'SWPC_SWMF_052811_2/raw_output/',
        'SWPC_iono': storage + 'SWPC_SWMF_052811_2/IONO-2D_CDF/',

        'TESTANALYTIC_cdf' : storage + 'TESTANALYTIC/GM/IO2/',
        'TESTANALYTIC_derived' : storage + 'TESTANALYTIC-derived/',

        'IMP10_RUN_SAMPLE_cdf' : storage + 'IMP10_RUN_SAMPLE/GM/IO2/',
        'IMP10_RUN_SAMPLE_derived' : storage + 'IMP10_RUN_SAMPLE-derived/',

        'interpolator': base + 'interpolators/'
        }

if conf['interpolator'] not in sys.path:
    sys.path.append(conf['interpolator'])

if base + 'util/' not in sys.path:
    sys.path.append(base + 'util/')

if base + 'physics/' not in sys.path:
    sys.path.append(base + 'physics/')

if base + 'derivatives/' not in sys.path:
    sys.path.append(base + 'derivatives/')

#if base + 'swmf/' not in sys.path:
#    sys.path.append(base + 'swmf/')

if base + 'timeseries/' not in sys.path:
    sys.path.append(base + 'timeseries/')

if base + 'cutplane/' not in sys.path:
    sys.path.append(base + 'cutplane/')

if conf['interpolator'] + 'kameleon/lib/python2.7/site-packages/ccmc/' not in sys.path:
    sys.path.append(conf['interpolator'] + 'kameleon/lib/python2.7/site-packages/ccmc/')

#if base + 'regions/' not in sys.path:
#    sys.path.append(base + 'regions/')

#if base + 'vtk/' not in sys.path:
#    sys.path.append(base + 'vtk/')

#if kameleon not in sys.path:
#    sys.path.append(kameleon)

#if kameleon + 'ccmc/' not in sys.path:
#    sys.path.append(kameleon + 'ccmc/')

#if not os.path.exists(base +'.logs/'):
#    os.makedirs(base +'.logs/')
#    print('created directory ' + base +'.logs/')

#if not os.path.exists(conf['DIPTSUR2_derived']+'timeseries/slices/'):
#    os.makedirs(conf['DIPTSUR2_derived']+'timeseries/slices/')
#    print('created directory '+conf['DIPTSUR2_derived']+'timeseries/slices/')

#if not os.path.exists(conf['SWPC_cdf']):
#    os.makedirs(conf['SWPC_cdf'])
#    print('Created directory ' + conf['SWPC_cdf'])
#
#if not os.path.exists(conf['SWPC_derived']):
#    os.makedirs(conf['SWPC_derived'])
#    print('Created directory ' + conf['SWPC_derived'])

#if not os.path.exists(conf['SCARR5_cdf']):
#    os.makedirs(conf['SCARR5_cdf'])
#    print('Created directory ' + conf['SCARR5_cdf'])

#if not os.path.exists(conf['SCARR5_derived']):
#    os.makedirs(conf['SCARR5_derived'])
#    print('Created directory ' + conf['SCARR5_derived'])

#if not os.path.exists(conf['SCARR1_cdf']):
#    os.makedirs(conf['SCARR1_cdf'])
#    print('Created directory ' + conf['SCARR1_cdf'])

#if not os.path.exists(conf['SCARR1_derived']):
#    os.makedirs(conf['SCARR1_derived'])
#    print('Created directory ' + conf['SCARR1_derived'])

#if not os.path.exists(conf['CARR_IMPULSE_cdf']):
#    os.makedirs(conf['CARR_IMPULSE_cdf'])
#    print('Created directory ' + conf['CARR_IMPULSE_cdf'])

#if not os.path.exists(conf['CARR_IMPULSE_derived']):
#    os.makedirs(conf['CARR_IMPULSE_derived'])
#    print('Created directory ' + conf['CARR_IMPULSE_derived'])
