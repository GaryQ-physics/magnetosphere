# write_VTKs

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config_paths import config
conf = config()

import B_field_lines_write
import structured_grid_write
import cut_plane
import longitude_lines_write

line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
time = line_list[0:5] + [0.]
Event = time + line_list[5:7]
T = tuple(time)
filename = conf["run_path"] + '3d__var_3_e' \
            + '%04d%02d%02d-%02d%02d%02d-000' % T + '.out.cdf'

cut_plane.writedata(Event)
B_field_lines_write.writevtk(Event)
longitude_lines_write.writevtk(Event)


if False:
    structured_grid_write.writevtk(Event, 'p')
    structured_grid_write.writevtk(Event, 'jy')
    structured_grid_write.writevtk(Event, 'dB')
    structured_grid_write.writevtk(Event, 'dBx')
    structured_grid_write.writevtk(Event, 'dBy')
    structured_grid_write.writevtk(Event, 'p')
    structured_grid_write.writevtk(Event, 'dBy', calcTotal=True)

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
            cut_plane.writedata(Event)
            B_field_lines_write.writevtk(Event)
            structured_grid_write.writevtk(Event, 'p')
            structured_grid_write.writevtk(Event, 'jy')
            structured_grid_write.writevtk(Event, 'dB_dV')
            structured_grid_write.writevtk(Event, 'dBlon_dV')
            structured_grid_write.writevtk(Event, 'dBlat_dV')
        else:
            print('ERROR: ' + filename + ' does not exitst')
    
    g.close()
