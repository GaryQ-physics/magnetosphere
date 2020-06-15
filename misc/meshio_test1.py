"""
pip install meshio

only runs with python3
"""

# meshio_test1

import numpy as np
import meshio
#import pathlib

in_fname = '/home/gary/magnetosphere/test_data/mesh_import.vtk'
out_fname = '/home/gary/magnetosphere/test_data/mesh_export.vtk'

mesh1 = meshio.read(
    in_fname,  # string, os.PathLike, or a buffer/open file
    file_format="vtk"  # optional if filename is a path; inferred from extension
)
"""
when using python 2, gives error: Unknown file format 'vtk' 
    if the optional file_format isn't given, gives error: Only VTK UNSTRUCTURED_GRID supported (not STRUCTURED_GRID)
"""

# mesh.points, mesh.cells, mesh.cells_dict, ...

print(mesh1.points)
print(mesh1.cells)

# mesh.vtk.read(filename) # is also possible


#meshio.write(out_fname, mesh1, file_format='vtk') # **kwargs

if False:
    points = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        ])
    cells = [
        ("triangle", np.array([[0, 1, 2]]))
    ]
    '''
    meshio.write_points_cells(
        out_fname,
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
        )
    '''
    
    mesh2 = meshio.Mesh(points, cells)
    meshio.write(
        out_fname,  # str, os.PathLike, or buffer/ open file
        mesh2,
        file_format="vtk",  # optional if first argument is a path; inferred from extension
    )
    print(mesh2.points)

if True:
    m = meshio.Mesh.read(in_fname, "vtk")  # same arguments as meshio.read
    m.write(out_fname, file_format="vtk")  # same arguments as meshio.write, besides `mesh`
