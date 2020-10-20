#deltaB_test.py

import os
os.system('f2py -c deltaB_subroutine.f90 -m deltaB_subroutine -DF2PY_REPORT_ON_ARRAY_COPY=1')

import numpy as np

import deltaB_subroutine

i=1

p1,p2 = deltaB_subroutine.rout(i)

ret = deltaB_subroutine.rout(i)

print(i)
print(ret)


case = 4
print('testing case %d'%(case))

if case==1:
    a = 2.*np.ones((7, 3), order='C')
    b = 3.*np.ones((7, 3), order='C')
    ax = a.__array_interface__['data'][0]

    deltaB_subroutine.arrtest(a, b, 7, 3)
    c = deltaB_subroutine.arrtest(a, b)

    print(a.__array_interface__['data'][0] == ax)
    print(c.__array_interface__['data'][0] == ax)
    print(c)
    print(a)

if case==2:
    a = 2.*np.ones((7, 3), order='F')
    b = 3.*np.ones((7, 3), order='F')
    ax = a.__array_interface__['data'][0]

    deltaB_subroutine.arrtest(a, b, 7, 3)
    c = deltaB_subroutine.arrtest(a, b)

    print(a.__array_interface__['data'][0] == ax)
    print(c.__array_interface__['data'][0] == ax)
    print(c)
    print(a)

if case==3:
    a = 2.*np.ones((7, 3), order='C')
    b = 3.*np.ones((7, 3), order='C')
    ax = a.__array_interface__['data'][0]

    a = a.transpose()
    b = b.transpose()

    print(np.isfortran(a))
    print(a.__array_interface__['data'][0] == ax)

    deltaB_subroutine.arrtest(a, b, 3, 7)
    c = deltaB_subroutine.arrtest(a, b)

    a = a.transpose()
    b = b.transpose()
    c = c.transpose()


    print(not np.isfortran(a))
    print(a.__array_interface__['data'][0] == ax)
    print(c.__array_interface__['data'][0] == ax)
    print(c)
    print(a)

if case==4:
    a = 2.*np.ones((4,5,3))
    b = 3.*np.ones((4,5,3))
    ax = a.__array_interface__['data'][0]
    print(not np.isfortran(a))

    a = a.transpose()
    b = b.transpose()
    deltaB_subroutine.cross(a, b, 5, 4)
    a = a.transpose()
    b = b.transpose()

    print(not np.isfortran(a))
    print(a.__array_interface__['data'][0] == ax)
    print(a)
    print(b)

