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
    soft, hard = int(4*2**30), int(4*2**30)
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
    if '-s' in sys.argv:
        storage = '/home/gary/media_sunspot/'

elif os.path.exists('/home/gquaresi/'):
    base = '/home/gquaresi/magnetosphere/'
    storage = '/media/sunspot/git-data/sblake/'
    #https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
    soft, hard = 72*2**30, 72*2**30
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
else:
    assert(False)


conf = {
        'storage': storage,
        'run_url': 'http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/',
        'run_path': storage + 'data/SCARR5_GM_IO2/IO2/',
        'run_path_derived': storage + 'SCARR5_GM_IO2-derived/',
        'base': base,
        'SWPC_raw': storage + 'SWPC_SWMF_052811_2/raw_output/',
        'SWPC_cdf_path': storage + 'SWPC_SWMF_052811_2/GM_CDF/',

        'SWPC_cdf': storage + 'SWPC_SWMF_052811_2/GM_CDF/', # = SWPC_cdf_path
        'SWPC_derived': storage + 'SWPC_SWMF_052811_2-derived/',

        'SCARR5_cdf': storage + 'SCARR5/GM/IO2/', # = run_path
        'SCARR5_derived': storage + 'SCARR5-derived/', # = run_path_derived
        'SCARR5_magfile': storage + 'SCARR5/MAG_FILES/',

        'SCARR1_cdf': storage + 'SCARR1/tmp/',
        'SCARR1_magfile': storage + 'SCARR1/MAG_FILES/',
        'SCARR1_derived': storage + 'SCARR1-derived/',

        'CARR_IMPULSE_cdf': storage + 'CARR_IMPULSE_3D/GM/IO2/',
        'CARR_IMPULSE_derived': storage + 'CARR_IMPULSE_3D-derived/',

        'DIPTSUR2_cdf' : storage + 'DIPTSUR2/GM/IO2/',
        'DIPTSUR2_derived' : storage + 'DIPTSUR2-derived/',

        'LUHMANN1979_cdf' : storage + 'LUHMANN1979/',
        'LUHMANN1979_derived' : storage + 'LUHMANN1979-derived/',

        'mag_server_url': 'http://mag.gmu.edu/git-data/sblake/',

        'SWPC_raw': storage + 'SWPC_SWMF_052811_2/raw_output/',
        'SWPC_iono': storage + 'SWPC_SWMF_052811_2/IONO-2D_CDF/',

        #'SCARR5_raw': storage + 'SCARR5_GM_IO2/IO2/',
        #'SCARR5_iono': storage + 'SCARR5_GM_IO2/IO2/',

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

if base + 'vtk/' not in sys.path:
    sys.path.append(base + 'vtk/')

if base + 'cutplane/' not in sys.path:
    sys.path.append(base + 'cutplane/')

#if kameleon not in sys.path:
#    sys.path.append(kameleon)

#if kameleon + 'ccmc/' not in sys.path:
#    sys.path.append(kameleon + 'ccmc/')

#if base + 'misc/' not in sys.path:
#    sys.path.append(base + 'misc/')

if base + 'regions/' not in sys.path:
    sys.path.append(base + 'regions/')

if not os.path.exists(conf['SWPC_cdf']):
    os.makedirs(conf['SWPC_cdf'])
    print('Created directory ' + conf['SWPC_cdf'])

if not os.path.exists(conf['SWPC_derived']):
    os.makedirs(conf['SWPC_derived'])
    print('Created directory ' + conf['SWPC_derived'])

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
