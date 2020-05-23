# magnetic_local_time (MLT)

#directory where magnetosphere folder is
m_path = '/home/gary/'


import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
#sys.path.append(conf["m_path"] + 'magnetosphere/data/')
import pos_sun as ps

path = conf["m_path"] + 'magnetosphere/data/'
fname = 'LOCALIZED_MLT.txt'
N=50

# units
deg = (np.pi/180.)
amin = deg/60.
hr = 1.
minn = hr/60.
s = minn/60.

f = open(path+fname,'w')
g = open(path+'LOCALIZED.txt','r')
print(g.readline()) # prints header line
line=0
for n in range(N):
    line=line+1
    # get string of the nth line of data
    line_string = g.readline()
    # get list of numbers from string  
    line_list = [float(i) for i in line_string.split()] # line_list "=" [year,month,day,hours,minutes,mlon,mlat]

    # get data
    Time = [line_list[0], line_list[1], line_list[2], line_list[3], line_list[4], 0.]
    MLONdeg = line_list[5]
    MLATdeg=line_list[6]
    phi=MLONdeg*deg
    #theta=np.pi/2. - MLAT*deg

    subsol_pt = ps.GSMtoMAG([1,0,0],Time,'car','car')
    phi0 = np.arctan2(subsol_pt[1], subsol_pt[0])
    
    delta = phi-phi0
    if(delta > np.pi):
        delta=delta-2.*np.pi
    elif(delta <= -np.pi):
        delta=delta+2.*np.pi
    MLT = 12.*hr + delta*24.*hr/(360*deg)

    #x=np.cos(phi)*np.sin(theta)
    #y=np.sin(phi)*np.sin(theta)
    #z=np.cos(theta)

    print "line ", line
    print(line_list)
    print(MLT)
    #print(subsol_pt)
    #print "GSM: ", ps.MAGtoGSM([x,y,z],month,day,year,UT)
    #print(delta)
    print Time
    print "Mdip", ps.MAGtoGSM([0.,0.,1.],Time,'car','car')
    # write to out file
    f.write(str(MLT/hr) + '\n')

g.close()
f.close()
