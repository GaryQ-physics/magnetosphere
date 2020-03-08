#!/home/gary/magnetosphere/kameleon/bin/python 

import sys
sys.path.append('../../../../lib/python2.7/site-packages/ccmc/')
import _CCMC as ccmc

kameleon = ccmc.Kameleon()
#kameleon.open('/home/gary/3d__var_3_e20031120-070000-000.out.cdf')
#interpolator = 0
#interpolator = kameleon.createNewInterpolator()
#print("interpreter created")
#kameleon.close()

def k_open(name):
	kameleon.open(name)
	print("file", name,"opened")
	#interpolator = kameleon.createNewInterpolator()

def main(argv):
	#print("main ran")
	if (len(argv) == 5):
	    filename = argv[0]
	    variable = argv[1]

	    c0 = float(argv[2])
	    c1 = float(argv[3])
	    c2 = float(argv[4])

	    #kameleon = ccmc.Kameleon()
	    #kameleon.open(filename)
	    kameleon.loadVariable(variable)
	    interpolator = kameleon.createNewInterpolator()
	    var = interpolator.interpolate(variable,c0, c1, c2)

	    # print variable, var

	    #kameleon.close()
	    return var
	else:
		print('Usage: <filename> <variable> x, y, z \n python kameleon_test rho -40 0 0')

def k_close():
	kameleon.close()
	print("file closed")

if __name__ == '__main__':
    main(sys.argv[1:])
    print("print when kp run as main")
else:
	print("Print when kp imported")
