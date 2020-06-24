# Requires Python 3.5+
import os
import sys
import imageio # Also need to install imageio-ffmpeg
from pathlib import Path

if not sys.version_info[0] == 3 and sys.version_info[1] >= 5:
    print("Python >= 3.5 required.")
    sys.exit(1)

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

vars = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']
delta = 0.1

for var in vars:
    image_dir = conf['run_path_derived'] + "cutplanes/" + var + "/"
    image_path = Path(image_dir)
    
    print('Reading list of files in ' + image_dir)
    images = list(image_path.glob('*' + '-' + var + '-' + '{0:.3f}'.format(delta) + '.png'))
    images.sort()
    
    print('Reading {0:d} files '.format(len(images)))
    image_list = []
    k = 0
    for file_name in images:
        if k > 100:
            break
        image_list.append(imageio.imread(file_name))
        fn = file_name.name
        print("{0} {1:s}".format(image_list[-1].shape, fn.replace(image_dir, "")))
        k = k + 1
    
    print('Generating mp4 using {0:d} files'.format(len(image_list)))
    outfile = conf['run_path_derived'] + "cutplanes/" + \
                    var + '-{0:.3f}'.format(delta) + '.mp4'
    print('Writing ' + outfile)
    imageio.mimwrite(outfile, image_list)
    print('Wrote ' + outfile)