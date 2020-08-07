
import numpy as np
import biot_savart_kameleon_interpolated_grid as bsk


global_xlims = (-224., 32.)
global_ylims = (-128., 128.)
global_zlims = (-128., 128.)


data = np.array([[2003, 11, 20, 7, 0, 57.50, 176.00]])

time = data[0, 0:5]
mlat = data[0, 5]
mlon = data[0, 6]
Event = data[0, :]

outs = np.nan*np.empty((6,))
for i in range(6):
    print('i = ' + str(i))
    if i == 0:
        outs[i] = 0.
        continue

    xlims = tuple(np.array(global_xlims)/float(2**i))
    ylims = tuple(np.array(global_ylims)/float(2**i))
    zlims = tuple(np.array(global_zlims)/float(2**i))

    outvect = bsk.integrate(time, mlat, mlon, para=True,
            Nx=None, xlims=xlims, dx=0.125,
            Ny=None, ylims=ylims, dy=0.125,
            Nz=None, zlims=zlims, dz=0.125,
            print_output=True)

    outs[i] = np.linalg.norm(outvect)

f = open('biot_savart_kameleon_demo2_outputfile.txt','w')
f.write('biot_savart_kameleon_demo2_outputfile\n')
np.savetxt(f, outs)
f.write('end\n')
f.close()

print(outs)
