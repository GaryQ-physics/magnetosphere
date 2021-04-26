#https://superuser.com/questions/592464/how-do-i-merge-a-folder-of-png-images-to-a-mp4-movie-using-avconv
#(python3.7) gary@gary-Inspiron-5567:~/media_sunspot/DIPTSUR2-derived/timeseries/slices$ ffmpeg -framerate 10 -pattern_type glob -i "xzplane_derivsb1_*.png" xzplane_ZZZ_ffconvert.mp4

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

def main(run):
    image_dir = '/media/sunspot/git-data/sblake/%s-derived/timeseries/slices/'%(run)
    image_path = Path(image_dir)
    
    print('Reading list of files in ' + image_dir)
    file_pattern = 'xzplane_derivsb1_*.png'
    image_files = list(image_path.glob(file_pattern))
    if len(image_files) == 0:
        print('No files found that match pattern ' + file_pattern)
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
    outfile = './xzplane_derivsb1.mp4'
    print('Writing ' + outfile)
    imageio.mimwrite(outfile, image_list)
    print('Wrote ' + outfile)

if __name__ == '__main__':
    main('DIPTSUR2')
