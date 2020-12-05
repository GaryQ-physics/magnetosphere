# adapted from https://gitlab.kitware.com/vtk/vtk/blob/32ae6c76b70c8da05cec8cabe3ccda8988388561/Filters/General/Testing/Python/streamTracer.py
# see ~/home/gary/VTKData/streamTracer.py
import os
import sys
import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from probe import probe
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

"""
trace()
	callable=None
		trace([[xmin, xmax, dx],[ymin, ymax, dy],[zmin, zmax, dz]], [bx,by,bz], points, method='scipy')
			Uses ??

		trace([[xmin, xmax, dx],[ymin, ymax, dy],[zmin, zmax, dz]], [bx,by,bz], points, method='vtk')
			Creates structured grid and passes to streamTracer

		trace([x,y,z],[bx,by,bz],points,method='scipy')
			Uses nninterpolator

		trace([x,y,z],[bx,by,bz],points,method='vtk')
			Use triangulate and then passes to streamTracer

	callable=function
		trace([[xmin, xmax, dx],[ymin, ymax, dy],[zmin, zmax, dz]], [bx,by,bz], points, method='scipy')
			Uses ??

		trace([[xmin, xmax, dx],[ymin, ymax, dy],[zmin, zmax, dz]], [bx,by,bz], points, method='vtk')
			Not supported

		trace([x,y,z],[bx,by,bz],points,method='scipy')
			Uses nninterpolator

		trace([x,y,z],[bx,by,bz],points,method='vtk')
			Not supported
	
traceFile (no callable option provided)
	traceFile('file.cdf', method='kameleon')
		Uses kamelon interpolator

	traceFile('file.cdf', method='scipy')
		Uses NNinterpolator callable

	traceFile('file.cdf', method='vtk')
		Uses triangulate to create grid

	traceFile('file.out', method='kameleon')
		Not supported

	traceFile('file.out', method='scipy')
		Uses NNinterpolator callable

	traceFile('file.out', method='vtk')
		Uses triangulate to create grid

"""
'''
trace(x,y,z,bx,by,bz,points,method='scipy',callable=None)
trace(x,y,z,bx,by,bz,points,method='vtk',callable=None)

(mostly done) trace(file.npy,method='scipy',callable=None)
1. trace(file.npy or file.vtk,method='vtk',callable=None) # Not ideal, but generally useful

(mostly done) trace(vtkObject,method='scipy',callable=vtkInterpolator)
(mostly done) trace(vtkObject,method='scipy',callable=SciPylnterpolator) (requires converting vtkObject to NumPy)
(mostly done) trace(vtkObject,method='vtk',callable=None)

2. trace(callable,method='scipy') # (Use for magnetic field models)
trace(callable,method='vtk') -> Won't work

(essentially done) trace(file.cdf,method='kameleon') # Uses native grid for interpolation. Fast.
(essentially done) trace(file.cdf,method='scipy',callable=probe) -> Slow b/c probe opens file each iteration

(essentially done) trace(file.cdf,method='scipy',callable=None) # reads on fine grid, uses SciPyInterpolator internally. 
(will be done when 1. is done) trace(file.cdf,method='vtk',callable=None) # reads on fine grid

Eventually want 
3. trace(file.out,method='vtk',callable=None) # reads on native grid, passes grid to vtk for tracing

May also want (for other models)
4. trace(file.cdf,method='vtk',callable=None) # reads on native grid, passes grid to vtk for tracing
'''


'''

trace(IC, Field, region=None, method='scipy')

-1. trace(IC, foo, region=None, method='vtk') -> ValueError
0. trace(IC, vtkObject, region=None, method='scipy') -> ValueError


1. trace(IC, foo, region=None, method='scipy') -> use scipy integrate on callable foo

2. trace(IC, array, region = {dictionary}, method='scipy') -> use scipy interpolator to make callable field, then use scipy integrate

3. trace(IC, file.cdf, region = {dictionary}, method='scipy') -> use probe on regular grid with scipy interpolator to make callable field, then use scipy integrate


4. trace(IC, vtk_object, region=None, method='vtk') -> uses vtk integrator on that vtk object

5. trace(IC, file.vtk, region=None, method='vtk') -> read vtk file to make vtk object, then use vtk integrator

6. trace(IC, file.cdf, region={dictionary}, method='vtk') -> use probe on regular grid to make vtk object, then use vtk integrator


7. trace(IC, file.cdf, region=None, method='kameleon') -> use kameleon integrator

(eventually want) 8. trace(IC, file.out, region=None, method='vtk') -> use native grid to make vtk object or vtk integrate

'''

def trace(IC, Field, region=None, method='scipy', debug=True):
    ''' Field = file.npy, file.vtk, file.cdf, tkObject, foo, 
       method = 'scipy','vtk'  '''
    import types
    sign = +1

    if method=='scipy':
        from scipy.integrate import odeint
        from scipy.interpolate import RegularGridInterpolator

        if isinstance(Field, types.FunctionType):
            def dXds(X, s):
                F = Field(X)
                Fmag = np.linalg.norm(F)
                if 1e-9 < Fmag < 1e+7:
                    return (sign/Fmag)*F
                return [0., 0., 0.]

            s_grid = np.arange(0., 10., 0.1)
            max_iterations = 100

            if IC.shape == (3,):
                IC = [IC]
            ret = []
            linenum = 0
            for X0 in list(IC):
                if debug:
                    print('linenum = ' + str(linenum))
                done = False
                solns = np.empty((0, 3)) # Combined solutions
                i = 0
                while not done:
                    if debug:
                        print('i = ' + str(i))
                    soln = odeint(dXds, X0, s_grid)
                    R = soln[:, 0]**2+soln[:, 1]**2 + soln[:, 2]**2
                    # define condition on the field line points
                    # Find first location where soln steps out-of-bounds
                    #tr = np.where( False == (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.) )        
                    # Boolean array.


                    tr = (R >= 1) & (soln[:,0] > -30.) & (np.abs(soln[:, 2]) < 20.)
                    # RuntimeWarning: invalid value encountered in greater_equal


                    # Indices where stop conditions satisfied
                    tr_out = np.where(tr == False)
                    if debug:
                        print(tr)
                    if tr_out[0].size > 0:
                        # Stop condition found at least once. Use solution up to that point.s
                        solns = np.vstack((solns, soln[0:tr_out[0][0] + 1, :]))
                        done = True
                    elif max_iterations == i + 1:
                        solns = np.vstack((solns, soln))   # return soln   faster?
                        done = True
                    else:
                        # New initial condition is stop point.
                        X0 = soln[-1, :]
                        # Append solution but exclude last value, which is the
                        # new initial condition.
                        solns = np.vstack((solns, soln[0:-1, :]))
                    i = i + 1
                ret.append(solns)
                linenum += 1

            if len(ret) == 1:
                return ret[0]
            return ret

        elif isinstance(Field, np.ndarray):
            assert(region is not None)

            from make_grid import make_axes
            ax_list = make_axes(region['xlims'], region['ylims'], region['zlims'], region['d'])

            # https://stackoverflow.com/questions/21836067/interpolate-3d-volume-with-numpy-and-or-scipy
            Fx_interp = RegularGridInterpolator(tuple(ax_list), Field[0, :,:,:])
            Fy_interp = RegularGridInterpolator(tuple(ax_list), Field[1, :,:,:])
            Fz_interp = RegularGridInterpolator(tuple(ax_list), Field[2, :,:,:])

            def Fcallable(v):
                return np.array([Fx_interp(v)[0], Fy_interp(v)[0], Fz_interp(v)[0]])

            return trace(IC, Fcallable, region=region, method=method)


        elif isinstance(Field, str):
            assert(region is not None)
            if Field[-4:] != '.cdf': raise ValueError
            fieldvar = 'b' # to generalize ?

            from make_grid import make_axes

            ax_list = make_axes(region['xlims'], region['ylims'], region['zlims'], region['d'])
            Nx = ax_list[0].size
            Ny = ax_list[1].size
            Nz = ax_list[2].size
            G2, G1, G3 = np.meshgrid(ax_list[1], ax_list[0], ax_list[2]) # different than in make_grid

            P = np.column_stack( (G1.flatten(), G2.flatten(), G3.flatten()) )
            Fx = probe(Field, P, var=fieldvar+'x', library='kameleon').reshape(Nx, Ny, Nz)
            Fy = probe(Field, P, var=fieldvar+'y', library='kameleon').reshape(Nx, Ny, Nz)
            Fz = probe(Field, P, var=fieldvar+'z', library='kameleon').reshape(Nx, Ny, Nz)

            return trace(IC, np.array([Fx, Fy, Fz]), region=region, method=method)

        else: raise ValueError

    elif method=='vtk':
        print(Field)
        print(type(Field))

        if isinstance(Field, vtk.vtkObject):

            rk = vtk.vtkRungeKutta45()
            # Create source for streamtubes
            streamer = vtk.vtkStreamTracer()
            streamer.SetInputConnection(Field.GetOutputPort())
            streamer.SetStartPosition(X[0],X[1],X[2])#SetSourceData()
            streamer.SetMaximumPropagation(20) ###
            #streamer.SetIntegrationStepUnit(2) # apperars overiden by next lines, see https://vtk.org/doc/nightly/html/classvtkStreamTracer.html#afe365e81e110f354065f5adc8401d589
            streamer.SetMinimumIntegrationStep(0.01)
            streamer.SetMaximumIntegrationStep(0.1)
            streamer.SetInitialIntegrationStep(0.01)
            #streamer.SetIntegrationDirection()
            #streamer.SetIntegrationDirectionToForward()
            #streamer.SetIntegrationDirectionToBackward()
            streamer.SetIntegrationDirectionToBoth()
            streamer.SetIntegrator(rk)
            streamer.SetRotationScale(0.5)
            streamer.SetMaximumError(1.0e-8)

            #https://stackoverflow.com/questions/38504907/reading-a-vtk-polydata-file-and-converting-it-into-numpy-array
            ## either order works ##
            polydata = streamer.GetOutput()
            streamer.Update()

            return dsa.WrapDataObject(polydata).Points

        elif isinstance(Field, str):
            if Field[-4:] == '.vtk':
                #make vtk object

                # read data
                reader = vtk.vtkStructuredGridReader()
                reader.SetFileName(Field)
                reader.Update()
                #force a read to occur

                return trace(IC, reader, region=region, method=method)

            elif Field[-4:] == '.cdf':
                assert(region is not None)
                # make vtk object
                return trace(IC, vtk_object, region=region, method=method)

            else: raise ValueError

        else: raise ValueError

    elif method=='kameleon':
        assert(Field[-4:] == '.cdf')
        pass #TODO later

    else: raise ValueError


def trace_old(X, view=False, write=False):
    datafile = "/tmp/tmp.vtk" # from ./magnetovis_vtk_demo.py

    # read data
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(datafile)
    reader.Update()
    #force a read to occur

    rk = vtk.vtkRungeKutta45()
    # Create source for streamtubes
    streamer = vtk.vtkStreamTracer()
    streamer.SetInputConnection(reader.GetOutputPort())
    streamer.SetStartPosition(X[0],X[1],X[2])#SetSourceData()
    streamer.SetMaximumPropagation(20) ###
    #streamer.SetIntegrationStepUnit(2) # apperars overiden by next lines, see https://vtk.org/doc/nightly/html/classvtkStreamTracer.html#afe365e81e110f354065f5adc8401d589
    streamer.SetMinimumIntegrationStep(0.01)
    streamer.SetMaximumIntegrationStep(0.1)
    streamer.SetInitialIntegrationStep(0.01)
    #streamer.SetIntegrationDirection()
    #streamer.SetIntegrationDirectionToForward()
    #streamer.SetIntegrationDirectionToBackward()
    streamer.SetIntegrationDirectionToBoth()
    streamer.SetIntegrator(rk)
    streamer.SetRotationScale(0.5)
    streamer.SetMaximumError(1.0e-8)

    print('#####################')
    print(streamer)
    print(type(streamer))
    print('#####################')
    print(streamer.GetOutputPort())
    print('#####################')

    if write:
        writer = vtk.vtkDataSetWriter()
        writer.SetInputConnection(streamer.GetOutputPort())
        writer.SetFileName("out.vtk")
        writer.Write()

    if view:
        outline = vtk.vtkStructuredGridOutlineFilter()
        outline.SetInputConnection(reader.GetOutputPort())
        mapOutline = vtk.vtkPolyDataMapper()
        mapOutline.SetInputConnection(outline.GetOutputPort())
        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(mapOutline)
        outlineActor.GetProperty().SetColor(0,0,0)

        mapStream = vtk.vtkPolyDataMapper()
        print(streamer.GetOutputPort())
        print(type(streamer.GetOutputPort()))
        mapStream.SetInputConnection(streamer.GetOutputPort())
        mapStream.SetScalarRange(reader.GetOutput().GetScalarRange())
        streamActor = vtk.vtkActor()
        streamActor.SetMapper(mapStream)

        return [outlineActor, streamActor, streamer]


    #https://stackoverflow.com/questions/38504907/reading-a-vtk-polydata-file-and-converting-it-into-numpy-array

    ## either order works ##
    polydata = streamer.GetOutput()
    streamer.Update()

    return dsa.WrapDataObject(polydata).Points



    #https://searchcode.com/file/25190165/Filters/Core/Testing/Python/fieldToRGrid.py/
    #ds2do = vtk.vtkDataSetToDataObjectFilter()
    #ds2do.SetInputConnection(streamer.GetOutputPort())
    #writer = vtk.vtkDataObjectWriter()
    #writer.SetInputConnection(ds2do.GetOutputPort())

    #ds_att = vtk.vtkDataSetAttributes()
    #ds_att.SetInputConnection(streamer.GetOutputPort())

    #pnts = vtk.vtkPoints()
    #pnts.SetData(streamer.GetOutputPort())
    #print(pnts)

    #rdr = vtk.vtkDataSetReader()
    #rdr.SetInputConnection(streamer.GetOutputPort())

def main():
    ren1 = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    # this test has some wireframe geometry
    # Make sure multisampling is disabled to avoid generating multiple
    # regression images
    # renWin SetMultiSamples 0
    renWin.AddRenderer(ren1)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    X = np.array([0.1,2.1,0.5])
    outlineActor, streamActor, streamer = trace_old(X, view=True)

    #print(streamer.GetOutputPort())

    ren1.AddActor(outlineActor)
    ren1.AddActor(streamActor)

    ren1.SetBackground(0.4,0.4,0.5)
    cam = ren1.GetActiveCamera()
    cam.SetPosition(20.,20.,20.) ###
    cam.SetFocalPoint(0.,0.,0.) ###
    cam.SetViewUp(0.311311,0.279912,0.908149)
    cam.SetClippingRange(1.12294,16.6226)
    renWin.SetSize(900,600)
    iren.Initialize()

    ### added ###
    renWin.Render()
    iren.Start()

    #print(streamer.GetOutputPort())

#vtkAlgorithmOutput

if __name__=='__main__':
    X = np.array([0.1,2.1,0.5])
    if False:
        main()
        arr = trace_old(X, view=False, write=True)
        print(arr)

    def foo(v):
        return np.array([1.,1.,1.])

    pm = 31.875
    reg =  {'xlims': (-pm, pm),
            'ylims': (-pm, pm),
            'zlims': (-pm, pm),
            'd': 0.25
            }
    cdffile = '/home/gary/magnetosphere/data/DIPTSUR2/GM/IO2/3d__var_2_e20190902-063000-000.out.cdf'

    #ret = trace(X, foo, region=None, method='scipy')
    #ret = trace(X, cdffile, region=reg, method='scipy')
    ret = trace(X, '/tmp/tmp.vtk', region=None, method='vtk')

    print('###################################')
    print(ret)

    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.plot3D(ret[:,0], ret[:,1], ret[:,2], '.')
    plt.show()
