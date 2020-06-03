# Requires Python 3.5+
import os
import sys
import imageio # Also need to install imageio-ffmpeg
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

v = 'p'
ext = '-' + v
image_path = Path(conf['run_path_derived'] + "cutplanes/")
images = list(image_path.glob('*' + ext + '.png'))
image_list = []
for file_name in images:
    image_list.append(imageio.imread(file_name))

print('Generating mp4 using {0:d} files'.format(len(image_list)))
outfile = conf['run_path_derived'] + "cutplanes/" + '3d__var_3_e' + v + '.mp4'
imageio.mimwrite(outfile, image_list)