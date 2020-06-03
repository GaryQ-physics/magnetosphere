# Generate cut plane for each file

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

debug = True

import cut_plane_plot as cp
from urlretrieve import urlretrieve

urlretrieve(conf['run_url'] + 'ls-1.txt', conf['run_path'] + 'ls-1.txt')

with open(conf['run_path'] + 'ls-1.txt','r') as f:
    files = f.readlines()

k = 0
times = []
for file in files:
    if k > 5:
        break
    k = k + 1
    filename = file.rstrip()
    if not filename[0:9] == '3d__var_3':
        continue

    tstr = filename[11:] 
    y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
    h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
    f = int(tstr[16:19])
    times.append([y, m, d, h, M, s, f])

    fname_full = conf['run_path'] + filename
    if not os.path.exists(fname_full):
        fileurl = conf['run_url'] + filename
        if debug:
            print('Downloading ' + fileurl)
            print('to')
        fname_tmp = fname_full + ".tmp"
        if debug:
            print(fname_tmp)
        urlretrieve(fileurl, fname_tmp)
        if debug:
            print('Downloaded ' + fileurl)
        os.rename(fname_tmp, fname_full)
        if debug:
            print('Renamed.')

    v = 'p'
    ext = '-' + v
    filename_png = conf["run_path_derived"] + "cutplanes/" + filename + ext + '.png'
        
    if os.path.exists(fname_full) and not os.path.exists(filename_png):
        cp.plot(times[-1], 'xz', None, v, 
                dx=0.5, dy=0.5,
                xlim=[-30, 30], ylim=[-20, 20], zlim=[0, 10],
                png=True, pngfile=filename_png, debug=debug)
