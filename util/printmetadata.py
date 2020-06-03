
import spacepy.pybats.bats as bats

############################################################################
# read in the 3d magnetosphere
d = '/Users/robertweigel/git/students/gquaresi/magnetosphere/data/SCARR5_GM_IO2/IO2/'

f = d + "3d__var_3_e20031120-070000-000.out"
print("Opening " + f)

data3d = bats.Bats2d(f)

# look at keys:
print(data3d.keys())

print(data3d.attrs)
print(data3d.meta)

print('------')

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import _CCMC as ccmc

d = '/Users/robertweigel/git/students/gquaresi/magnetosphere/data/SCARR5_GM_IO2/IO2/'
f = d + '3d__var_3_e20031120-070000-000.out.cdf'

kameleon = ccmc.Kameleon()
print("Opening " + f)
kameleon.open(f)

for i in range(kameleon.getNumberOfGlobalAttributes()):
    gname = kameleon.getGlobalAttributeName(i)
    gattr = kameleon.getGlobalAttribute(gname)
    if gname != 'README':
        print gname, gattr.toString()
        
for i in range(kameleon.getNumberOfVariables()):
    varname  = kameleon.getVariableName(i)
    min_attr = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
    max_attr = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
    units = kameleon.getVisUnit(varname)
    units2 = kameleon.getNativeUnit(varname)
    print varname, '\t', min_attr,'\t', max_attr, units, units2

kameleon.close()    