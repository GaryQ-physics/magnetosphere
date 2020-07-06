"""
This script (WILL) compare the dB value calculated by biot-savart integral in
this repository to that pulled from mag_grid_files
"""

import read_mag_grid_files as rmg

if False:
    rmg.analyzedata('mag_grid_e20031120-070000.out', 176.00, 57.50)
else:
    rmg.analyzedata((2003,11,20,7,0,0), 176.00, 57.50)
