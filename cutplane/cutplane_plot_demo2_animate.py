# Requires Python 3.5+, imageio, imageio-ffmpeg
import os
import sys
import imageio
from pathlib import Path

if not sys.version_info[0] == 3 and sys.version_info[1] >= 5:
    print("Python >= 3.5 required.")
    sys.exit(1)

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

run = 'SCARR5'

# TODO: Generate list based on subdirectories (and exclude directory minmax.)
variables = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']
#variables = ['p']
delta = 0.125
plot_type = 1
pattern = 'type_{0:d}_delta_{1:.3f}'.format(plot_type, delta)
#pattern = '{0:.3f}'.format(delta)

for var in variables:
    image_dir = conf[run + '_derived'] + "cutplanes/" + var + "/"
    image_path = Path(image_dir)
    
    print('Reading list of files in ' + image_dir)
    file_pattern = '*' + '-' + var + '-' + pattern + ".png"
    image_files = list(image_path.glob(file_pattern))
    if len(image_files) == 0:
        print('No files found that match pattern ' + file_pattern)
        continue
    image_files.sort()
    
    print('Reading {0:d} files '.format(len(image_files)))
    image_list = []
    k = 0
    for image_file in image_files:
        image_list.append(imageio.imread(image_file))
        file_name = image_file.name
        print("{0:d}/{1:d} size = {2}, {3:s}" \
              .format(k+1, len(image_files), 
                      image_list[-1].shape, 
                      file_name.replace(image_dir, "")))
        k = k + 1
    
    print('Generating mp4 using {0:d} files'.format(len(image_list)))
    outfile = conf[run + '_derived'] + "cutplanes/" + \
                    var + "-" + pattern + '.mp4'
    print('Writing ' + outfile)
    imageio.mimwrite(outfile, image_list)
    print('Wrote ' + outfile)
