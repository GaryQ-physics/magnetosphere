"""
pip install pyvista (only works python3)

Create an Unstructured Grid from NumPy arrays.
save as ./.created_grid.vtk and also plot it with pyvistas included plotter
"""

import pyvista as pv
import numpy as np

import vtk
HEXAHEDRON_int = vtk.VTK_HEXAHEDRON ## this is an integer whose value is vtk's convention
del vtk

# these points will all be shared between the cells
points = np.array([[0. , 0. , 0. ],
                   [1. , 0. , 0. ],
                   [0.5, 0. , 0. ],
                   [1. , 1. , 0. ],
                   [1. , 0.5, 0. ],
                   [0. , 1. , 0. ],
                   [0.5, 1. , 0. ],
                   [0. , 0.5, 0. ],
                   [0.5, 0.5, 0. ],
                   [1. , 0. , 0.5],
                   [1. , 0. , 1. ],
                   [0. , 0. , 0.5],
                   [0. , 0. , 1. ],
                   [0.5, 0. , 0.5],
                   [0.5, 0. , 1. ],
                   [1. , 1. , 0.5],
                   [1. , 1. , 1. ],
                   [1. , 0.5, 0.5],
                   [1. , 0.5, 1. ],
                   [0. , 1. , 0.5],
                   [0. , 1. , 1. ],
                   [0.5, 1. , 0.5],
                   [0.5, 1. , 1. ],
                   [0. , 0.5, 0.5],
                   [0. , 0.5, 1. ],
                   [0.5, 0.5, 0.5],
                   [0.5, 0.5, 1. ]])


# Each cell in the cell array needs to include the size of the cell
# and the points belonging to the cell.  In this example, there are 8
# hexahedral cells that have common points between them.
cells = np.array([[ 8,  0,  2,  8,  7, 11, 13, 25, 23],
                  [ 8,  2,  1,  4,  8, 13,  9, 17, 25],
                  [ 8,  7,  8,  6,  5, 23, 25, 21, 19],
                  [ 8,  8,  4,  3,  6, 25, 17, 15, 21],
                  [ 8, 11, 13, 25, 23, 12, 14, 26, 24],
                  [ 8, 13,  9, 17, 25, 14, 10, 18, 26],
                  [ 8, 23, 25, 21, 19, 24, 26, 22, 20],
                  [ 8, 25, 17, 15, 21, 26, 18, 16, 22]]).ravel()

# each cell is a VTK_HEXAHEDRON
celltypes = np.empty(8, dtype=np.uint8)
celltypes[:] = HEXAHEDRON_int

# note the `cells` array is 1 dimensional.
# the `offset` array will point to the start of each cell in the `cells` array
# the offset array is not needed for vtk DataFile Version 9.0 or newer
offset = np.array([ 0, 9, 18, 27, 36, 45, 54, 63])

# Effectively, when visualizing a VTK unstructured grid, it will
# sequentially access the cell array by first looking at each index of
# cell array (based on the offset array), and then read the number of
# points based on the first value of the cell.  In this case, the
# VTK_HEXAHEDRON is described by 8 points.

# for example, the 5th cell would be accessed by vtk with:
start_of_fifthcell = offset[4]
n_points_in_fifthcell = cells[start_of_fifthcell]
indices_in_fifthcell = cells[start_of_fifthcell + 1: start_of_fifthcell + n_points_in_fifthcell + 1]
print(indices_in_fifthcell)

if True:
    #if you are using or newer, you do not need to input the offset array:
    grid = pv.UnstructuredGrid(cells, celltypes, points)
elif True:
    # if you are not using VTK 9.0 or newer, you must use the offset array
    grid = pv.UnstructuredGrid(offset, cells, celltypes, points)
else:
    # Alternate versions:
    grid = pv.UnstructuredGrid({HEXAHEDRON_int: cells.reshape([-1, 9])[:, 1:]}, points)
    grid = pv.UnstructuredGrid({HEXAHEDRON_int: np.delete(cells, np.arange(0, cells.size, 9))}, points)

grid.save('./.created_grid.vtk', binary=False)
# plot the grid (and suppress the camera position output)
_ = grid.plot(show_edges=True)
