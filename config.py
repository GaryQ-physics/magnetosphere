import os
import sys

if os.path.exists('/Users/robertweigel/'):
    base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'
    kameleon = '/Users/robertweigel/kameleon/lib/python2.7/site-packages/'
if os.path.exists('/Users/weigel/'):
    base = '/Users/weigel/git/magnetosphere/'
    kameleon = '/Users/weigel/kameleon/lib/python2.7/site-packages/'
elif os.path.exists('/home/weigel/') and False:
    base = '/home/weigel/git/magnetosphere/'
    kameleon = '/home/weigel/kameleon/lib/python2.7/site-packages/'
elif os.path.exists('/home/gary/'):
    base = '/home/gary/magnetosphere/'
    kameleon = '/home/gary/magnetosphere/kameleon/lib/python2.7/site-packages/'
    storage = base
elif os.path.exists('/home/gquaresi/'):
    base = '/home/gquaresi/magnetosphere/'
    kameleon = '/home/gquaresi/'
    storage = '/media/solar-backup/tmp/'
else:
    assert(False)


conf = {
        'data_path': base + 'data/',
        'run_url': 'http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/',
        'run_path': storage + 'data/SCARR5_GM_IO2/IO2/',
        'run_path_derived': storage + 'data/SCARR5_GM_IO2-derived/',
        'base': base,
        'SWPC_raw': storage + 'data/SWPC_SWMF_052811_2/raw_output/',
        'SWPC_cdf_path': storage + 'data/SWPC_SWMF_052811_2/GM_CDF/',

        'SWPC_cdf': storage + 'data/SWPC_SWMF_052811_2/GM_CDF/', # = SWPC_cdf_path
        'SWPC_derived': storage + 'data/SWPC_SWMF_052811_2-derived/',
        'SCARR5_cdf': storage + 'data/SCARR5_GM_IO2/IO2/', # = run_path
        'SCARR5_derived': storage + 'data/SCARR5_GM_IO2-derived/', # = run_path_derived
        'SCARR1_magfile': storage + 'data/SCARR5_GM_IO2/IO2/',
        'SCARR1_cdf': storage + 'data/SCARR1/tmp/',
        'SCARR1_magfile': storage + 'data/SCARR1/MAG_FILES/',
        'SCARR1_derived': storage + 'data/SCARR1-derived/',
        'mag_server_url': 'http://mag.gmu.edu/git-data/sblake/',

        'SWPC_raw': storage + 'data/SWPC_SWMF_052811_2/raw_output/',
        'SWPC_iono': storage + 'data/SWPC_SWMF_052811_2/IONO-2D_CDF/',

        'SCARR5_raw': storage + 'data/SCARR5_GM_IO2/IO2/',
        'SCARR5_iono': storage + 'data/SCARR5_GM_IO2/IO2/'
    }

if base + 'util/' not in sys.path:
    sys.path.append(base + 'util/')

if base + 'physics/' not in sys.path:
    sys.path.append(base + 'physics/')

if base + 'vtk/' not in sys.path:
    sys.path.append(base + 'vtk/')

if base + 'cutplane/' not in sys.path:
    sys.path.append(base + 'cutplane/')

if kameleon not in sys.path:
    sys.path.append(kameleon)

if kameleon + 'ccmc/' not in sys.path:
    sys.path.append(kameleon + 'ccmc/')

if base + 'misc/' not in sys.path:
    sys.path.append(base + 'misc/')

if not os.path.exists(conf['SWPC_cdf']):
    os.makedirs(conf['SWPC_cdf'])
    print('Created directory ' + conf['SWPC_cdf'])

if not os.path.exists(conf['SWPC_derived']):
    os.makedirs(conf['SWPC_derived'])
    print('Created directory ' + conf['SWPC_derived'])

if not os.path.exists(conf['SCARR5_cdf']):
    os.makedirs(conf['SCARR5_cdf'])
    print('Created directory ' + conf['SCARR5_cdf'])

if not os.path.exists(conf['SCARR5_derived']):
    os.makedirs(conf['SCARR5_derived'])
    print('Created directory ' + conf['SCARR5_derived'])

if not os.path.exists(conf['SCARR1_cdf']):
    os.makedirs(conf['SCARR1_cdf'])
    print('Created directory ' + conf['SCARR1_cdf'])

if not os.path.exists(conf['SCARR1_derived']):
    os.makedirs(conf['SCARR1_derived'])
    print('Created directory ' + conf['SCARR1_derived'])
