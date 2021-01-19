import os
import sys
import numpy as np

from config import conf

sys.path.append(conf['base'] + 'kameleonV-0.2.3/')
import kameleonV


#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00002000_n0002973.out.cdf'


#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_130700_000.cdf' # 405.5 kB
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_060700_000.cdf' # 2.7 MB
filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf' # 2.7 MB

arr = np.array([10.,10.])
#P = np.column_stack([arr,arr,arr])

out = kameleonV.interpolate(filename, arr, arr, arr, 'y')
print(out)
