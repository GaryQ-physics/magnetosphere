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
    if grid == 'POLYDATA'
        num_points = points.shape[0]
        f.write('\nPOINTS '+str(num_points)+' float\n')

    if ftype=='BINARY':
        points = np.array(points, dtype='>f')
        f.write(points.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, points)

    f.write('\n')

    if grid == 'POLYDATA':
        if connectivity == 'triangulation':
            assert(ftype == 'ASCII')
            extra_circle = 1
            closing_points = 2
            head_count = 1
            extra_point = 1
            points_on_circle = 0
            x_first_circle = X[0]
            for x_search in X:
                if x_search == x_first_circle:
                    points_on_circle += 1
                else:
                    break
            
            num_circles = int(num_points/points_on_circle)
            
            triangle_line = ''
            lines_of_triangles = num_circles - extra_circle
            head_triangle = 2 * points_on_circle + closing_points
            points_on_triangle_strips = (2 
                                        * points_on_circle 
                                        + closing_points 
                                        + head_count) * lines_of_triangles
            f.write('\n')
            f.write('TRIANGLE_STRIPS {} {}\n'.format(lines_of_triangles, 
                                                     points_on_triangle_strips))

            for i in range(num_points - points_on_circle):
                triangle_line = triangle_line + ' {} {}'.format(i, i+points_on_circle)
                

                if (i+1) % int(points_on_circle) == 0:
                    if i == points_on_circle: 
                        f_repeat_point = 0
                        l_repeat_point = f_repeat_point + points_on_circle
                    else: 
                        f_repeat_point = i - points_on_circle + extra_point
                        l_repeat_point = f_repeat_point + points_on_circle
                    line = '{} {} {} {}\n'.format(head_triangle, 
                                               triangle_line,
                                               f_repeat_point,
                                               l_repeat_point)
                    f.write(line)
                    triangle_line = ''
            f.write('\n')

    if grid == 'STRUCTURED_GRID':
        f.write('POINT_DATA ' + str(Nx*Ny*Nz) + '\n')
    if grid == 'POLYDATA':
        f.write('POINT_DATA ' + str(num_points) + '\n')

    if texture == 'SCALARS':
        f.write('SCALARS ' + point_data_name + ' float 1\n') # number with float???
        f.write('LOOKUP_TABLE default\n')
    if texture == 'VECTORS':
        f.write('VECTORS ' + point_data_name + ' float\n')
    if texture == 'TEXTURE_COORDINATES':
        f.write('TEXTURE_COORDINATES ' + point_data_name + ' 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping


    if ftype=='BINARY':
        point_data = np.array(point_data, dtype='>f')
        f.write(point_data.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, point_data)

    f.close()
    if debug:
        print("Wrote " + out_filename)
