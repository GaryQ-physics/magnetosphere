# config_paths
import os
def config():
    if os.path.exists('/Users/robertweigel/'):
        base = '/Users/robertweigel/git/students/gquaresi/'
        return {'k_path': '/Users/robertweigel/',
                'f_path': base + 'magnetosphere/data/',
                'm_path': base}
    else:
        return {'k_path': '/home/gary/magnetosphere/',
                'f_path': '/home/gary/magnetosphere/data/',
                'm_path': '/home/gary/'}
