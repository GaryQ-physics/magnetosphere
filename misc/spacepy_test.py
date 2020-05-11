# spacepy_test

import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config_paths import config
conf = config()
#import spacepy as sp
import spacepy.coordinates as sc
from spacepy.time import Ticktock
#spacepy.coordinates.Coords


cvals = sc.Coords([[1,20,40],[1,20,20]], 'MAG', 'sph')
print(cvals)
#print cvals.x
cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
newcoord = cvals.convert('GSM', 'car')
print newcoord
print newcoord.x

'''
cvals = sc.Coords(np.array([1,20,40]), 'MAG', 'sph')
print(cvals)
#print cvals.x
cvals.ticks = Ticktock('2002-02-02T12:00:00', 'ISO') # add ticks
newcoord = cvals.convert('GSM', 'car')
print newcoord
print newcoord.x[0]
'''
