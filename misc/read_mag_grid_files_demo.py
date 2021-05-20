"""
This script (WILL) compare the dB value calculated by biot-savart integral in
this repository to that pulled from mag_grid_files
"""

import read_mag_grid_files as rmg

if False:
    rmg.analyzedata('mag_grid_e20031120-070000.out', 57.50, 176.00)  #gives assertion error
elif False:
    rmg.analyzedata('/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/mag_grid_e20031120-070000.out', 57.50, 176.00)
    rmg.analyzedata('/home/gary/magnetosphere/data/SCARR1/MAG_FILES/mag_grid_e20031120-070000.out', 57.50, 176.00)
    ret = rmg.analyzedata((2003,11,20,7,0), 57.50, 176.00)
    print(ret)
else:
    ret = rmg.analyzedata('/home/gary/temp/mag_grid_e20190902-063000.out', 57.50, 176.00)
    print(ret)
    ret = rmg.analyzedata('/home/gary/temp/mag_grid_e20190902-063000.out', 5.7327586207E+01, 176.00)
    print(ret)
