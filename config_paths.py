# config_paths
import os
def config():
    #k_path = 'kameleon/lib/python2.7/site-packages/'
    if os.path.exists('/Users/robertweigel/'):
        return {'k_path': '/Users/robertweigel/', 'f_path': '/Users/robertweigel/Desktop', 'm_path': '/home/gary/'}
    else:
        return {'k_path': '/home/gary/magnetosphere/', 'f_path': '/home/gary/magnetosphere/events/', 'm_path': '/home/gary/'}
