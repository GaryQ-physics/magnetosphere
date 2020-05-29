# config_paths
import os
import sys

def config():
    if os.path.exists('/Users/robertweigel/'):
        base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'
        kameleon = '/Users/robertweigel/'
    else:
        base = '/home/gary/magnetosphere/'
        kameleon = base

    ret = {
            'data_path': base + 'data/',
            'run_path': base + 'data/SCARR5_GM_IO2/IO2/',
            'run_path_derived': base + 'data/SCARR5_GM_IO2-derived/'
        }
    sys.path.append(base + 'events/')
    sys.path.append(base + 'vtk/')
    sys.path.append(kameleon + 'kameleon/lib/python2.7/site-packages/')
    sys.path.append(kameleon + 'kameleon/lib/python2.7/site-packages/ccmc/')
    
    if not os.path.exists(ret['run_path_derived']):
        os.mkdir(ret['run_path_derived'])
        
    return ret
