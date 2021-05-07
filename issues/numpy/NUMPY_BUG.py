import numpy as np
import sys
print(sys.version)
print('numpy version : ' +str(np.__version__))
print('\n')

################ using python floats
gridspacing = 0.5
xlims = (0.0, 1.0)
ylims = (0.0, 1.0)
zlims = (0.0, 1.0)

grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                ylims[0]:ylims[1]+gridspacing:gridspacing,
                zlims[0]:zlims[1]+gridspacing:gridspacing ]

print('using python built in float for limits and spacing yields array of float:')
print(grid)

################# using numpy float64
gridspacing = np.float64(0.5)
xlims = (np.float64(0.0), np.float64(1.0))
ylims = (np.float64(0.0), np.float64(1.0))
zlims = (np.float64(0.0), np.float64(1.0))

grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                ylims[0]:ylims[1]+gridspacing:gridspacing,
                zlims[0]:zlims[1]+gridspacing:gridspacing ]

print('using numpy.float64 for limits and spacing yields array of floats')
print('that agrees with the python one:')
print(grid)

################# using numpy float32
gridspacing = np.float32(0.5)
xlims = (np.float32(0.0), np.float32(1.0))
ylims = (np.float32(0.0), np.float32(1.0))
zlims = (np.float32(0.0), np.float32(1.0))

grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                ylims[0]:ylims[1]+gridspacing:gridspacing,
                zlims[0]:zlims[1]+gridspacing:gridspacing ]

print('using numpy.float32 for limits and spacing yields array of integers')
print('that differs in value due to rounding:')
print(grid)
