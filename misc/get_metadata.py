# get_metadata.py
# python2

import sys
import os
import numpy as np

import spacepy.pybats.bats as bt

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import _CCMC as ccmc
from urlretrieve import urlretrieve

# read in the 3d magnetosphere
filename = "3d__var_3_e20031120-070000-000.out"
fname = conf['run_path'] + filename

# download if not found
if not os.path.exists(fname):
    fileurl = conf['run_url'] + filename
    print('Downloading ' + fileurl)
    fname_tmp = fname + ".tmp"
    urlretrieve(fileurl, fname_tmp)
    os.rename(fname_tmp, fname)
    print('Downloaded ' + fname)


# use pybats
data3d = bt.Bats2d(fname)
#print(data3d.keys()) #.keys() even though data3d not dictionary
for key in data3d.meta.keys():
    print key, '->', data3d.meta[key]
print('------------------')
for key in data3d.attrs.keys():
    print key, '->', data3d.attrs[key]

# use kameleon
kameleon = ccmc.Kameleon()
kameleon.open(fname + '.cdf')
print("Opened " + fname + '.cdf')
interpolator = kameleon.createNewInterpolator()


print(kameleon.getNumberOfGlobalAttributes())
print(kameleon.getNumberOfVariableAttributes())

NumAttributes = kameleon.getNumberOfGlobalAttributes() + kameleon.getNumberOfVariableAttributes()
print('------------------')
for i in range(NumAttributes):
    gname = kameleon.getGlobalAttributeName(i)
    gattr = kameleon.getGlobalAttribute(gname)
    if gname != 'README' and gattr.toString() != 'STRING: : ':
        print(i)
        print gname, gattr.toString()

print('------------------')   
for i in range(NumAttributes):
    vname = kameleon.getVariableAttributeName(i)
    vattr = kameleon.getVariableAttribute('jy',vname)
    if vattr.toString() != 'STRING: : ':
        print(i)
        print vname, vattr.toString()


kameleon.close()
print("Closed " + fname + '.cdf')
