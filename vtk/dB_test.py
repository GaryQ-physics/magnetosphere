# dB_test

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()

import structured_grid_write as sgw

line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
time = line_list[0:5] + [0, 0.]
Event = time + line_list[5:7]
T = tuple(time)
filename = conf["run_path"] + '3d__var_3_e' \
            + '%04d%02d%02d-%02d%02d%02d-%03d' % T + '.out.cdf'


vals=[]
mults=[]
for i in range(4):
    mult = (i+1)/2.
    val = sgw.Compute(Event, 'dB_EW', calcTotal=True, retTotal=True, mult = mult)
    vals.append(val)
    mults.append(mult)




import matplotlib.pyplot as plt
plt.plot(mults,vals)
plt.ylabel('dB_EW/nT')
plt.xlabel('mult (dx=0.3R_e/mult)')
plt.show()