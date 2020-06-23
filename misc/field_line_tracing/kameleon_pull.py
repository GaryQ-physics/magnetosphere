# uncompyle6 version 3.7.1
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.18 |Anaconda, Inc.| (default, Apr 23 2020, 22:42:48) 
# [GCC 7.3.0]
# Embedded file name: /home/gary/magnetosphere/kameleon/bin/ccmc/examples/python/kameleon_pull.py
# Compiled at: 2020-03-08 17:11:16
import sys
import os
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../../' )
from config import conf
import _CCMC as ccmc
kameleon = ccmc.Kameleon()

def k_open(name):
    kameleon.open(name)
    print ('file', name, 'opened')


def main(argv):
    if len(argv) == 5:
        filename = argv[0]
        variable = argv[1]
        c0 = float(argv[2])
        c1 = float(argv[3])
        c2 = float(argv[4])
        kameleon.loadVariable(variable)
        interpolator = kameleon.createNewInterpolator()
        var = interpolator.interpolate(variable, c0, c1, c2)
        return var
    print 'Usage: <filename> <variable> x, y, z \n python kameleon_test rho -40 0 0'


def k_close():
    kameleon.close()
    print 'file closed'


if __name__ == '__main__':
    main(sys.argv[1:])
    print 'print when kp run as main'
else:
    print 'Print when kp imported'
# okay decompiling kameleon_pull.pyc
