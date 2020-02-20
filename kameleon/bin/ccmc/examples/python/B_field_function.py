# B_field_function

#    ../../../python2 kameleon_test.py ~/magnetosphere/kameleon-plus-compiled/3d__var_3_e20031120-070000-000.out.cdf bx 5 0 0

import sys
#sys.path.append("/home/gary/kameleon/bin/ccmc/examples/python/")
import kameleon_pull as kp

#exec(open("~/kameleon/bin/ccmc/examples/python/kameleon_test.py").read())

#from importlib.machinery import SourceFileLoader

#foo = SourceFileLoader("kameleon_test.py", "/home/gary/kameleon/bin/ccmc/examples/python/").load_module()
#foo.MyClass()



import numpy as np

filename='/home/gary/magnetosphere/3d__var_3_e20031120-070000-000.out.cdf'

def Bx(x,y,z):
	ret=kp.main([filename,'bx',x,y,z])
	return ret

#end
if __name__ == '__main__':
    print("Executed when bfield invoked directly")
    print(Bx(1.,2.,3.))
else: 
    print("Executed when bfield imported")
