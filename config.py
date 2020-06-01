import os
import sys

if os.path.exists('/Users/robertweigel/'):
    base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'
    kameleon = '/Users/robertweigel/kameleon/lib/python2.7/site-packages/'
else:
    base = '/home/gary/magnetosphere/'
    kameleon = '/home/gary/magnetosphere/kameleon/lib/python2.7/site-packages/'

conf = {
        'data_path': base + 'data/',
        'run_path': base + 'data/SCARR5_GM_IO2/IO2/',
        'run_path_derived': base + 'data/SCARR5_GM_IO2-derived/'
    }

if base + 'events/' not in sys.path:
    sys.path.append(base + 'events/')

if base + 'vtk/' not in sys.path:
    sys.path.append(base + 'vtk/')
    
if kameleon not in sys.path:
    sys.path.append(kameleon)

if kameleon + 'ccmc/' not in sys.path:
    sys.path.append(kameleon + 'ccmc/')
    
if not os.path.exists(conf['run_path_derived']):
    os.mkdirs(conf['run_path_derived'])
    print('Created directory ' + conf['run_path_derived'])
    
