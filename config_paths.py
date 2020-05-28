# config_paths
import os
import sys

def config():
    if os.path.exists('/Users/robertweigel/'):
        kameleon = '/Users/robertweigel/'
        base = '/Users/robertweigel/git/students/gquaresi/'
        ret = {'run_path': base + 'magnetosphere/data/SCARR5_GM_IO2/IO2/',
               'run_path_derived': base + 'magnetosphere/data/SCARR5_GM_IO2-derived/'}
    else:
        base = '/home/gary/magnetosphere/'
        kameleon = base
        ret = {'run_path': base + 'data/SCARR5_GM_IO2/',
               'run_path_derived': base + 'data/SCARR5_GM_IO2-derived/'}

    sys.path.append(base + 'events/')
    sys.path.append(base + 'vtk/')
    sys.path.append(kameleon + 'kameleon/lib/python2.7/site-packages/')
    sys.path.append(kameleon + 'kameleon/lib/python2.7/site-packages/ccmc/')
    
    return ret
