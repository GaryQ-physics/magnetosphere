"""
This script (WILL) compare the dB value calculated by biot-savart integral in
this repository to that pulled from mag_grid_files
"""

import read_mag_grid_files as rmg

if False:
    rmg.analyzedata('mag_grid_e20031120-070000.out', 57.50, 176.00)
else:
    ret = rmg.analyzedata((2003,11,20,7,0), 57.50, 176.00)
    print(ret)