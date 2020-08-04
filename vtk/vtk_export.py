import numpy as np


def writevtk(out_filename, points, point_data, connectivity, texture,
             point_data_name = 'point data', title='Title', ftype='BINARY', grid='STRUCTURED_GRID', debug=True):
    """
    ftype = 'BINARY' or 'ASCII'
    grid = 'STRUCTURED_GRID' or 'TRIANGLE_STRIPS' or
    """

    if grid == 'STRUCTURED_GRID':
        Nx,Ny,Nz = connectivity

    assert(out_filename[0] == '/')
    f = open(out_filename,'w')
    if debug:
        print("Writing " + out_filename)
    f.write('# vtk DataFile Version 3.0\n')
    f.write(title + '\n')
    f.write(ftype + '\n')
    f.write('DATASET ' + grid + '\n')

    if grid == 'STRUCTURED_GRID':
        f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
        f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')

    if ftype=='BINARY':
        points = np.array(points, dtype='>f')
        f.write(points.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, points)

    f.write('\n')
    if grid == 'STRUCTURED_GRID':
        f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')

    if texture == 'SCALARS':
        f.write('SCALARS ' + point_data_name + ' float 1\n') # number with float???
        f.write('LOOKUP_TABLE default\n')
    if texture == 'VECTORS':
        f.write('VECTORS ' + point_data_name + ' float\n')
    if texture == 'TEXTURE_COORDINATES':
        f.write('TEXTURE_COORDINATES ' + point_data_name + ' 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping


    # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping

    '''
    f.write('TEXTURE_COORDINATES TextureCoordinates 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping


    if var=='J':
        f.write('VECTORS ' + var + ' float\n')
    else:
        f.write('SCALARS ' + var + ' float 1\n')
        f.write('LOOKUP_TABLE default\n')
    '''

    if ftype=='BINARY':
        point_data = np.array(point_data, dtype='>f')
        f.write(point_data.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, point_data)

    f.close()
    if debug:
        print("Wrote " + out_filename)
