# config_paths
import os
def config():
    if os.path.exists('/Users/robertweigel/'):
        base = '/Users/robertweigel/git/students/gquaresi/'
        return {'k_path': '/Users/robertweigel/',
                'f_path': base + 'magnetosphere/data/',
                'm_path': base,
                'run_path': base + 'magnetosphere/data/SCARR5_GM_IO2/IO2/',
                'run_path_derived': base + 'magnetosphere/data/SCARR5_GM_IO2-derived/'}
    else:
        base = '/home/gary/'
        return {'k_path': base + 'magnetosphere/',
                'f_path': base + 'magnetosphere/data/',
                'm_path': base,
                'run_path': base + 'magnetosphere/data/SCARR5_GM_IO2/',
                'run_path_derived': base + 'magnetosphere/data/SCARR5_GM_IO2-derived/'}
