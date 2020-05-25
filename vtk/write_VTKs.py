# write_VTKs

import sys
import os
import numpy as np
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
import field_line_vtk_script
import kameleon_structured_grid_write
import cut_plane

#MLONdeg = 68.50
#MLATdeg = 50.00
#Event = [2003, 11, 20, 7, 0, 0, 68.50, 50.00]

Event = [2003, 11, 20, 7, 0, 0, 241.00, 55.00] # [year,month,day,hours,minutes,seconds,MLONdeg,MLATdeg]

cut_plane.writevtk(Event)
field_line_vtk_script.writevtk(Event)
kameleon_structured_grid_write.writevtk(Event,'p')
kameleon_structured_grid_write.writevtk(Event,'jy')
kameleon_structured_grid_write.writevtk(Event,'dB_dV')
kameleon_structured_grid_write.writevtk(Event,'dBlon_dV')
kameleon_structured_grid_write.writevtk(Event,'dBlat_dV')
