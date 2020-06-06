import numpy as np
from niceticks import niceticks

ymin = -3.1
ymax = 1.1
n = 5
ticks = niceticks(-3.1, 1.1, n, debug=False)
print('ymin       = {0:f}'.format(ymin))
print('tickmin    = {0:f}'.format(ticks[0]))
print('ymax       = {0:f}'.format(ymax))
print('tickmax    = {0:f}'.format(ticks[-1]))
print('n          = {0:d}'.format(n))
print('len(ticks) = {0:d}'.format(len(ticks)))
print('ticks      = ' + ", ".join(np.array_str(ticks, precision=2, suppress_small=True).split()))
