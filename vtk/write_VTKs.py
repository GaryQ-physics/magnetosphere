# write_VTKs

import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

print 'c0'
#import b_field_lines_write
print 'c1'
import structured_grid_write
print 'c2'
import cutplane
print 'c3'
#import longitude_lines_write
#import j_field_lines_write
import earth_write
from events import eventlist


'''
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

#cut_plane.writedata(Event)
#cutplane.writedata(time, mlat, mlon)

#B_field_lines_write.writevtk(Event)
#longitude_lines_write.writevtk(Event) #need ad
#J_vector_field_write.writevtk(Event)
#earth_write.writevtk(Event) #need ad
#structured_grid_write.writevtk(Event, 'dB', binary=True, dx=0.2, dy=0.2, dz=0.2)
#structured_grid_write.writevtk(Event, 'jy', dx=0.2, dy=0.2, dz=0.2)
#structured_grid_write.writevtk(Event, 'dB_EW', calcTotal=True)
'''

line_list = [2003, 11, 20, 7, 0, 176.00, 57.50]

data = np.array([[2003, 11, 20, 7, 0, 57.50, ]])
#data = eventlist()

print('data.shape = ' + str(data.shape))

import time as tm
to = tm.time()

#time = np.array([2003, 11, 20, 7, 0])
#mlat = 57.50
#mlon = 176.00

time = np.array([2006, 12, 15, 7, 7, 0])
mlat = 57.50
mlon = 176.00

structured_grid_write.cdf_to_structured_grid('SWPC', time, mlat, mlon, 'p',
            xlims=(-56., 8.), ylims=(-32., 32.), zlims=(-32., 32.),
            d=0.25)

'''
N = 1
for i in range(N):
    time = data[i, 0:5]
    mlat = data[i, 5]
    mlon = data[i, 6]
    Event = data[i, :]
    print time, mlat, mlon


    #structured_grid_write.writevtk(Event, 'p')
    #longitude_lines_write.writevtk(Event)
    #earth_write.writevtk(Event)
    #b_field_lines_write.writevtk(Event)
    #j_field_lines_write.writevtk(Event)
'''



tf = tm.time()
print('run time = ' + str(tf-to))
#run time = 45.9153749943  with usekV=True
#run time = 60.5226249695  with usekV=False


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
