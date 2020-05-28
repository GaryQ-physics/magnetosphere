# write_VTKs

import os
import sys

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()

import field_line_vtk_script
import kameleon_structured_grid_write
import cut_plane


line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
time = line_list[0:5] + [0.]
Event = time + line_list[5:7]
T = tuple(time)
filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'
cut_plane.writevtk(Event)
field_line_vtk_script.writevtk(Event)
kameleon_structured_grid_write.writevtk(Event, 'p')
kameleon_structured_grid_write.writevtk(Event, 'jy')
kameleon_structured_grid_write.writevtk(Event, 'dB_dV')
kameleon_structured_grid_write.writevtk(Event, 'dBlon_dV')
kameleon_structured_grid_write.writevtk(Event, 'dBlat_dV')

if False:
    N = 10

    g = open(conf["run_path_derived"] + 'LOCALIZED.txt', 'r')
    print(g.readline()) # prints header line
    line = 0
    for n in range(N):
        line = line + 1
        # get string of the nth line of data
        line_string = g.readline()
        # get list of numbers from string  
        line_list = [float(i) for i in line_string.split()] # line_list "=" [year,month,day,hours,minutes,mlon,mlat]
    
        # get data
        time = line_list[0:5] + [0.]
        Event = time + line_list[5:7]
        T = tuple(time)
        filename = conf["run_path"] + '3d__var_3_e' + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'
        if os.path.exists(filename):
            cut_plane.writevtk(Event)
            field_line_vtk_script.writevtk(Event)
            kameleon_structured_grid_write.writevtk(Event, 'p')
            kameleon_structured_grid_write.writevtk(Event, 'jy')
            kameleon_structured_grid_write.writevtk(Event, 'dB_dV')
            kameleon_structured_grid_write.writevtk(Event, 'dBlon_dV')
            kameleon_structured_grid_write.writevtk(Event, 'dBlat_dV')
        else:
            print('ERROR: ' + filename + ' does not exitst')
    
    g.close()
