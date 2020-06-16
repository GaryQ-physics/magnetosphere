"""
pip install uvw

only runs with python3, but gives nonsense as is
"""

# uvw_test1

import numpy as np
from uvw import RectilinearGrid, DataArray

fname = '/home/gary/magnetosphere/test_data/uvw_test.vtk'
a=5
b=5
c=10

# Creating coordinates
x = np.linspace(-0.5, 0.5, a)
y = np.linspace(-0.5, 0.5, b)
z = np.linspace(-0.9, 0.9, c)

# Creating the file (with possible data compression)
grid = RectilinearGrid(fname, (x, y, z))

# A centered ball
x, y, z = np.meshgrid(x, y, z, indexing='ij')
r = np.sqrt(x**2 + y**2 + z**2)
ball = r < 0.3

# Some multi-component multi-dimensional data
data = np.zeros([a, b, c, 3, 3])
data[ball, ...] = np.array([[0, 1, 0],
                            [1, 0, 0],
                            [0, 1, 1]])

# Some cell data
cell_data = np.zeros([9, 19, 29])
cell_data[0::2, 0::2, 0::2] = 1

# Adding the point data (see help(DataArray) for more info)
grid.addPointData(DataArray(data, range(3), 'ball'))
# Adding the cell data
grid.addCellData(DataArray(cell_data, range(3), 'checkers'))
grid.write()
