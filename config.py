import os
import sys

if os.path.exists('/Users/robertweigel/'):
    base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'
    kameleon = '/Users/robertweigel/kameleon/lib/python2.7/site-packages/'
if os.path.exists('/Users/weigel/'):
    base = '/Users/weigel/git/magnetosphere/'
    kameleon = '/Users/weigel/kameleon/lib/python2.7/site-packages/'
elif os.path.exists('/home/weigel/'):
    base = '/home/weigel/git/magnetosphere/'
    kameleon = '/home/weigel/kameleon/lib/python2.7/site-packages/'
else:
    base = '/home/gary/magnetosphere/'
    kameleon = '/home/gary/magnetosphere/kameleon/lib/python2.7/site-packages/'

conf = {
        'data_path': base + 'data/',
        'run_url': 'http://mag.gmu.edu/git-data/sblake/SCARR5_GM_IO2/IO2/',
        'run_path': base + 'data/SCARR5_GM_IO2/IO2/',
        'run_path_derived': base + 'data/SCARR5_GM_IO2-derived/',
        'base': base
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

if not os.path.exists(conf['run_path_derived']):
    os.makedirs(conf['run_path_derived'])
    print('Created directory ' + conf['run_path_derived'])
