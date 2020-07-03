# write_VTKs

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

import B_field_lines_write
import J_field_lines_write
import structured_grid_write
import cut_plane
import cutplane
import longitude_lines_write
import earth_write
from events import events

line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]
mlat = line_list[6]
mlon = line_list[5]
time_orig = line_list[0:5]

time = time_orig + [0, 0.]
Event = time + line_list[5:7]
T = tuple(time)
filename = conf["run_path"] + '3d__var_3_e' \
            + '%04d%02d%02d-%02d%02d%02d-%03d' % T + '.out.cdf'

#structured_grid_write.writevtk(Event, 'J', calcTotal=True)
#ret = structured_grid_write.Compute(Event, 'dB_EW', calcTotal=True)
#structured_grid_write.writevtk(Event, 'dB_EW', binary=True, calcTotal=True)

#structured_grid_write.writevtk(Event, 'p', dx=1., dy=1., dz=1., fname='/tmp/testvtk.vtk')

#J_field_lines_write.writevtk(Event)

cut_plane.writedata(Event)
#cutplane.writedata(time, mlat, mlon)

#B_field_lines_write.writevtk(Event)
#longitude_lines_write.writevtk(Event) #need ad
#J_vector_field_write.writevtk(Event)
#earth_write.writevtk(Event) #need ad
#structured_grid_write.writevtk(Event, 'dB', binary=True, dx=0.2, dy=0.2, dz=0.2)
#structured_grid_write.writevtk(Event, 'jy', dx=0.2, dy=0.2, dz=0.2)
#structured_grid_write.writevtk(Event, 'dB_EW', calcTotal=True)

data = events()

N = 1
for i in range(N):
    time = data[i,0:5]
    mlat = data[5]
    mlon = data[6]
    Event = data[i,:]




if False:
    #J_field_lines_write.writevtk(Event)
    #cut_plane.writedata(Event)
    #B_field_lines_write.writevtk(Event)
    #longitude_lines_write.writevtk(Event) #need ad
    #J_vector_field_write.writevtk(Event)
    #earth_write.writevtk(Event) #need ad
    #structured_grid_write.writevtk(Event, 'jy')
    #structured_grid_write.writevtk(Event, 'dB_EW', calcTotal=True)
    structured_grid_write.writevtk(Event, 'dB', calcTotal=True)
    structured_grid_write.writevtk(Event, 'p')
    structured_grid_write.writevtk(Event, 'dB')
    structured_grid_write.writevtk(Event, 'dB_EW')
    structured_grid_write.writevtk(Event, 'dB_NS', calcTotal=True)

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
            structured_grid_write.writevtk(Event, 'dB')
            structured_grid_write.writevtk(Event, 'dB_EW')
            structured_grid_write.writevtk(Event, 'dB_NS')
        else:
            print('ERROR: ' + filename + ' does not exitst')
    
    g.close()
