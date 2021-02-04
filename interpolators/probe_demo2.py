# depreciated
import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from probe import probe

#filename ='/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-170100-085.out.cdf'
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00002000_n0002973.out.cdf'

#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_130700_000.cdf' # 405.5 kB
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/IONO-2D.swmf.it061215_060700_000.cdf' # 2.7 MB
#filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf' # 2.7 MB

if False:
    filename ='/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-170100-085.out.cdf'
    X = np.array([[-5, 0, 0], [-10, 0, 0]])

    J = probe(filename, X, var=['jx','jy','jz'], debug=False, dictionary=False, library='kameleonV')
    #J = probe(filename, X, var=['jx','jy','jz'], debug=False, dictionary=False, library='kameleon')

else:
    filename = '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf' # 2.7 MB
    X = np.array([[-5, 0], [-10, 0]])
    J = probe(filename, X, var=['jx','jy','jz'], debug=False, dictionary=False, library='pycdf')

print(J)
