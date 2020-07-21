import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import biot_savart_kameleon as bsk
import util
import cxtransform as cx
import read_ccmc_datafiles as r_ccmc

global_xlims = (-224., 32.)
global_ylims = (-128., 128.)
global_zlims = (-128., 128.)

data, headers = r_ccmc.getdata(2006)
times = np.array(data[:, 0:7], dtype=int)
index = 1100 # choose

time = times[index, :]
print('time = ' + str(time))
filename = util.time2SWPCfile(time)

outs = np.nan*np.empty((6,))
for i in range(6):
    print('i = ' + str(i))
    if i == 0 or i ==1 or i ==2 or i ==3 :
        outs[i] = 0.
        continue

    xlims = tuple(np.array(global_xlims)/float(2**i))
    ylims = tuple(np.array(global_ylims)/float(2**i))
    zlims = tuple(np.array(global_zlims)/float(2**i))

    YKClat = 62.480
    YKClon = 245.518

    mpos = cx.GEOtoMAG([1., YKClat, YKClon] , time, 'sph', 'sph')
    mlat = mpos[1]
    mlon = mpos[2]

    print(time)
    print(filename)

    outvect = bsk.run(time, mlat, mlon, para=True, filename=filename,
            Nx=None, xlims=xlims, dx=0.125,
            Ny=None, ylims=ylims, dy=0.125,
            Nz=None, zlims=zlims, dz=0.125,
            print_output=True)

    outs[i] = np.linalg.norm(outvect)

f = open('biot_savart_SWPC_outputfile.txt','w')
f.write('biot_savart_SWPC_outputfile\n')
f.write('time = ' + str(time) + '\n')
f.write('from kameleon:\n')
np.savetxt(f, outs)
f.write('from SWMF:\n')
f.write(str(headers[13:16]) + '\n')
f.write(str(data[index,13:16]) + '\n')
f.write('dB_norm = ' + str( np.linalg.norm(data[index,13:16]) ) + '\n')
f.write('end\n')
f.close()

print(outs)
print(str( np.linalg.norm(data[index,13:16]) ) )
