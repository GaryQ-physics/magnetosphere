import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
from printmetadata import printmetadata

printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out.cdf')
printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out')
