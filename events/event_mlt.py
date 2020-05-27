# Add MLT column to LOCALIZED.txt

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
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
    
    if line[0] == "#": # Header line
        newlines.append(line.rstrip() + "   MLT\n")
        continue    
    
    line_list = [float(i) for i in line.split()]

    # get data
    time = [line_list[0], line_list[1], line_list[2], line_list[3], line_list[4], 0.]
    MLONdeg = line_list[5]
    MLATdeg = line_list[6]

    MLT = ps.MLTfromMAG(MLONdeg, time)

    newlines.append(line.rstrip() + " {0:.2f}\n".format(MLT))

f = open(outfile, 'w')
print('Writing ' + outfile)
f.writelines(newlines)
print('Wrote ' + outfile)
