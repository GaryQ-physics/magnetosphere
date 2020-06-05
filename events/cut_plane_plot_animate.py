# Requires Python 3.5+
import os
import sys
import imageio # Also need to install imageio-ffmpeg
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

var = 'p'
delta = 0.5
image_path = Path(conf['run_path_derived'] + "cutplanes/")
images = list(image_path.glob('*' + '-' + var + '-' + '{0:.3f}'.format(delta) + '.png'))
images.sort()

image_list = []
for file_name in images:
    image_list.append(imageio.imread(file_name))

print('Generating mp4 using {0:d} files'.format(len(image_list)))
outfile = conf['run_path_derived'] + "cutplanes/" + \
            '3d__var_3_e' + '-' + '{0:.3f}'.format(delta) + '.mp4'
print('Writing ' + outfile)
imageio.mimwrite(outfile, image_list)
print('Wrote ' + outfile)