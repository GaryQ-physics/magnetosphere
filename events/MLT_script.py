# Add MLT column to LOCALIZED.txt

import sys
import os

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()

import pos_sun as ps

infile = conf["run_path_derived"] + 'LOCALIZED.txt'
outfile = conf["run_path_derived"] + 'LOCALIZED_MLT.txt'

print('Opening ' + infile)
g = open(infile,'r')
lines = g.readlines()
print('Opened ' + infile)
g.close()

print('Computing MLT for {0:d} events'.format(len(lines)))
newlines = []
for line in lines:    
    
    if line[0] == "#":
        newlines.append(line.rstrip() + "   MLT\n")
        continue # Skip header line
        
    # get list of numbers from string  
    line_list = [float(i) for i in line.split()]

    # get data
    time = [line_list[0], line_list[1], line_list[2], line_list[3], line_list[4], 0.]
    MLONdeg = line_list[5]
    MLATdeg = line_list[6]

    MLT = ps.MLTfromMAG(MLONdeg, time)

    #print(line_list)
    #print("MLT = {0:.2f}".format(MLT))

    newlines.append(line.rstrip() + " {0:,.2f}\n".format(MLT))
    #print(subsol_pt)
    #print "GSM: ", ps.MAGtoGSM([x,y,z],month,day,year,UT)
    #print(delta)
    #print "Mdip", ps.MAGtoGSM([0.,0.,1.],time,'car','car')
    # write to out file

#    f.write(str(MLT) + '\n')

f = open(outfile, 'w')
print('Writing ' + outfile)
f.writelines(newlines)
print('Wrote ' + outfile)
