import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cxtransform as cx

def events():
    """Event information
    
    Returns
    -------
    np.ndarray of [year, month, day, hour, minute, mag_lat, mag_lon, mlt]
    """
    infile = conf["run_path_derived"] + 'LOCALIZED.txt'

    data = np.genfromtxt(infile, skip_header=1)    

    mlt = cx.MAGtoMLT(data[:, 5], data[:, 0:5])

    # Swap mlat and mlon colums so in expected order (lat then long)
    data[:, [6,5]] = data[:, [5,6]]
    
    data = np.hstack((data, np.reshape(mlt, (mlt.shape[0], 1))))
    
    return data