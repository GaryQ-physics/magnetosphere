# config_paths
import os
def config():
    if os.path.exists('/Users/robertweigel/'):
        base = '/Users/robertweigel/git/students/gquaresi/magnetosphere/'
        return {'k_path': '/Users/robertweigel/',
                'f_path': base + '/data/',
                'm_path': base}
    else:
        return {'k_path': '/home/gary/magnetosphere/',
                'f_path': '/home/gary/magnetosphere/events/',
                'm_path': '/home/gary/'}
