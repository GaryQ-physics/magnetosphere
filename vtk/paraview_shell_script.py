# paraview_shell_script

###copy to shell
'''
execfile('paraview_shell_script.py')
'''
###

#fieldFile_path = m_path + 'magnetosphere/vtk/kameleon_field_paraview.vtk'

import sys
import os
import numpy as np
#print os.path.dirname(os.path.abspath(__file__)) + '/../'
#sys.path.append( os.path.dirname(os.path.abspath(_file_)) + '/../' )
sys.path.append( './../' )  # for some reason the other way wont work in paraview
from config import conf
#import pos_sun as ps
#import paraview.simple as pvs
from paraview.simple import *  #needed for paraview 5.8


Event = [2003, 11, 20, 7, 0, 0, 176.00, 57.50]
time = Event[0:6]
tag = '_%04d:%02d:%02dT%02d:%02d:%02d' % tuple(time)
subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])

Nlong=5
Nb = 6
N=Nb+1+Nlong
var='dB_EW'

cut_plane_name = 'cut_plane_info_%.2f_%.2f' %(Event[7], Event[6]) + tag + '.txt'
grid_name = 'structured_grid_' + var + tag + '.000' + '.vtk'
earth_name = 'earth' + tag + '.vtk'


#-----------------
#file = conf["data_path"] + "topography/world.topo.200401.3x5400x2700.png-ParaView.png"
file = conf["data_path"] + "topography/world.topo.2004{0:02d}.3x5400x2700.png".format(Event[1])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Legacy VTK Reader'
rotated_spherevtk = LegacyVTKReader(FileNames=[conf["run_path_derived"] + subdir + earth_name])
# show data in view
rotated_spherevtkDisplay = Show(rotated_spherevtk, renderView1)
# trace defaults for the display properties.
rotated_spherevtkDisplay.Representation = 'Surface'
rotated_spherevtkDisplay.ColorArrayName = [None, '']
rotated_spherevtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
rotated_spherevtkDisplay.SelectOrientationVectors = 'None'
rotated_spherevtkDisplay.ScaleFactor = 0.4
rotated_spherevtkDisplay.SelectScaleArray = 'None'
rotated_spherevtkDisplay.GlyphType = 'Arrow'
rotated_spherevtkDisplay.GlyphTableIndexArray = 'None'
rotated_spherevtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
rotated_spherevtkDisplay.PolarAxes = 'PolarAxesRepresentation'
rotated_spherevtkDisplay.ScalarOpacityUnitDistance = 0.15493986305312726

texProxy = servermanager.CreateProxy("textures", "ImageTexture")
texProxy.GetProperty("FileName").SetElement(0, file)
texProxy.UpdateVTKObjects()
rotated_spherevtkDisplay.Texture = texProxy

Render()
#--------------


f = open(conf["run_path_derived"] + subdir + cut_plane_name,'r')

# get string of the 1st line of data (Mdipole components)
Mdipole_string = f.readline()
# get list of components from string  
Mdipole = [float(i) for i in Mdipole_string.split()]
# make array
Mdipole=np.array(Mdipole)

#next line is U1, ect
U1_string = f.readline()
U1 = [float(i) for i in U1_string.split()]
U1=np.array(U1)
U2_string = f.readline()
U2 = [float(i) for i in U2_string.split()]
U2=np.array(U2)
U3_string = f.readline()
U3 = [float(i) for i in U3_string.split()]
U3=np.array(U3)

# create a new 'Legacy VTK Reader'
structured_gridvtk = LegacyVTKReader(FileNames=[conf["run_path_derived"] + subdir + grid_name])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [932, 802]

# show data in view
structured_gridvtkDisplay = Show(structured_gridvtk, renderView1)
# trace defaults for the display properties.
structured_gridvtkDisplay.Representation = 'Outline'
structured_gridvtkDisplay.ColorArrayName = ['POINTS', '']
structured_gridvtkDisplay.OSPRayScaleArray = var
structured_gridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
structured_gridvtkDisplay.SelectOrientationVectors = 'None'
structured_gridvtkDisplay.ScaleFactor = 21.0
structured_gridvtkDisplay.SelectScaleArray = var
structured_gridvtkDisplay.GlyphType = 'Arrow'
structured_gridvtkDisplay.GlyphTableIndexArray = var
structured_gridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
structured_gridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
structured_gridvtkDisplay.ScalarOpacityUnitDistance = 5.766431907004487

# reset view to fit data
renderView1.ResetCamera()

# change representation type
structured_gridvtkDisplay.SetRepresentationType('Points')

# create a new 'Slice'
slice1 = Slice(Input=structured_gridvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-95.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# get color transfer function/color map for var
point_scalarsLUT = GetColorTransferFunction(var)

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', var]
slice1Display.LookupTable = point_scalarsLUT
slice1Display.OSPRayScaleArray = var
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 2.0
slice1Display.SelectScaleArray = var
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = var
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)


# create a new 'Slice'
slice2 = Slice(Input=structured_gridvtk)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [-95.0, 0.0, 0.0]
slice2.SliceType.Normal = [0.,0.,1.]

# get color transfer function/color map for var
point_scalarsLUT = GetColorTransferFunction(var)

# show data in view
slice2Display = Show(slice2, renderView1)
# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['POINTS', var]
slice2Display.LookupTable = point_scalarsLUT
slice2Display.OSPRayScaleArray = var
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 2.0
slice2Display.SelectScaleArray = var
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = var
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice2Display.SetScalarBarVisibility(renderView1, True)



# create a new 'Slice'
slice3 = Slice(Input=structured_gridvtk)
slice3.SliceType = 'Plane'
slice3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [0.0, 0.0, 0.0]
slice3.SliceType.Normal = [1.,0.,0.]

# get color transfer function/color map for var
point_scalarsLUT = GetColorTransferFunction(var)

# show data in view
slice3Display = Show(slice3, renderView1)
# trace defaults for the display properties.
slice3Display.Representation = 'Surface'
slice3Display.ColorArrayName = ['POINTS', var]
slice3Display.LookupTable = point_scalarsLUT
slice3Display.OSPRayScaleArray = var
slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice3Display.SelectOrientationVectors = 'None'
slice3Display.ScaleFactor = 2.0
slice3Display.SelectScaleArray = var
slice3Display.GlyphType = 'Arrow'
slice3Display.GlyphTableIndexArray = var
slice3Display.DataAxesGrid = 'GridAxesRepresentation'
slice3Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice3Display.SetScalarBarVisibility(renderView1, True)



# create a new 'Slice'
slice4 = Slice(Input=structured_gridvtk)
slice4.SliceType = 'Plane'
slice4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Origin = [0.0, 0.0, 0.0]
slice4.SliceType.Normal = U3

# get color transfer function/color map for var
point_scalarsLUT = GetColorTransferFunction(var)

# show data in view
slice4Display = Show(slice4, renderView1)
# trace defaults for the display properties.
slice4Display.Representation = 'Surface'
slice4Display.ColorArrayName = ['POINTS', var]
slice4Display.LookupTable = point_scalarsLUT
slice4Display.OSPRayScaleArray = var
slice4Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice4Display.SelectOrientationVectors = 'None'
slice4Display.ScaleFactor = 2.0
slice4Display.SelectScaleArray = var
slice4Display.GlyphType = 'Arrow'
slice4Display.GlyphTableIndexArray = var
slice4Display.DataAxesGrid = 'GridAxesRepresentation'
slice4Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice4Display.SetScalarBarVisibility(renderView1, True)


# hide data in view
Hide(structured_gridvtk, renderView1)



h = 15.
#------------------
# x axis
#------------------
cylinderX = Cylinder()
cylinderX.Radius = 0.05
cylinderX.Center = [0.0, h/2., 0.0]
cylinderX.Height = h
cylinderXDisplay = Show(cylinderX, renderView1)
cylinderXDisplay.ColorArrayName = [None, '']
cylinderXDisplay.DiffuseColor = [1.0, 0.0, 0.0]
# Default Cylinder is orientated along y-axis.
# To get x cylinder, rotate by -90 around z-axis.
cylinderXDisplay.Orientation = [0.0, 0.0, -90.0]
# cone x
coneX = Cone()
# Properties modified on coneX
coneX.Resolution = 30
coneX.Radius = 0.2
coneX.Height = 0.4
coneX.Center = [h, 0.0, 0.0]
coneX.Direction = [1., 0., 0.]
# show data in view
coneXDisplay = Show(coneX, renderView1)
# trace defaults for the display properties.
coneXDisplay.Representation = 'Surface'
coneXDisplay.ColorArrayName = [None, '']
coneXDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coneXDisplay.SelectOrientationVectors = 'None'
coneXDisplay.ScaleFactor = 0.04000000357627869
coneXDisplay.SelectScaleArray = 'None'
coneXDisplay.GlyphType = 'Arrow'
coneXDisplay.GlyphTableIndexArray = 'None'
coneXDisplay.DataAxesGrid = 'GridAxesRepresentation'
coneXDisplay.PolarAxes = 'PolarAxesRepresentation'
# change solid color
coneXDisplay.DiffuseColor = [1.0, 0.0, 0.0]
for i in range(int(h)-1):
    # create a new 'Sphere'
    sph = Sphere()

    # Properties modified on sph
    sph.Center = [i+1, 0.0, 0.0]
    sph.Radius = 0.2
    sph.ThetaResolution = 10
    sph.PhiResolution = 10

    # show data in view
    sphDisplay = Show(sph, renderView1)
    # trace defaults for the display properties.
    sphDisplay.Representation = 'Surface'
    sphDisplay.ColorArrayName = [None, '']
    sphDisplay.OSPRayScaleArray = 'Normals'
    sphDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    sphDisplay.SelectOrientationVectors = 'None'
    sphDisplay.ScaleFactor = 0.2
    sphDisplay.SelectScaleArray = 'None'
    sphDisplay.GlyphType = 'Arrow'
    sphDisplay.GlyphTableIndexArray = 'None'
    sphDisplay.DataAxesGrid = 'GridAxesRepresentation'
    sphDisplay.PolarAxes = 'PolarAxesRepresentation'

    # change solid color
    sphDisplay.DiffuseColor = [1.0, 0.0, 0.0]








#------------------
# y axis
#------------------
cylinderY = Cylinder()
cylinderY.Radius = 0.05
cylinderY.Center = [0.0, h/2., 0.0]
cylinderY.Height = h
cylinderYDisplay = Show(cylinderY, renderView1)
cylinderYDisplay.ColorArrayName = [None, '']
cylinderYDisplay.DiffuseColor = [1.0, 1.0, 0.5]
# Default Cylinder is orientated along y-axis, so 
# the following statement is not needed.
cylinderYDisplay.Orientation = [0.0, 0.0, 0.0]
# cone y
coneY = Cone()
# Properties modified on coneY
coneY.Resolution = 30
coneY.Radius = 0.2
coneY.Height = 0.4
coneY.Center = [0.0, h, 0.0]
coneY.Direction = [0., 1., 0.]
# show data in view
coneYDisplay = Show(coneY, renderView1)
# trace defaults for the display properties.
coneYDisplay.Representation = 'Surface'
coneYDisplay.ColorArrayName = [None, '']
coneYDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coneYDisplay.SelectOrientationVectors = 'None'
coneYDisplay.ScaleFactor = 0.04000000357627869
coneYDisplay.SelectScaleArray = 'None'
coneYDisplay.GlyphType = 'Arrow'
coneYDisplay.GlyphTableIndexArray = 'None'
coneYDisplay.DataAxesGrid = 'GridAxesRepresentation'
coneYDisplay.PolarAxes = 'PolarAxesRepresentation'
# change solid color
coneYDisplay.DiffuseColor = [1.0, 1.0, 0.5]
for i in range(int(h)-1):
    # create a new 'Sphere'
    sph = Sphere()

    # Properties modified on sph
    sph.Center = [0., i+1, 0.]
    sph.Radius = 0.2
    sph.ThetaResolution = 10
    sph.PhiResolution = 10

    # show data in view
    sphDisplay = Show(sph, renderView1)
    # trace defaults for the display properties.
    sphDisplay.Representation = 'Surface'
    sphDisplay.ColorArrayName = [None, '']
    sphDisplay.OSPRayScaleArray = 'Normals'
    sphDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    sphDisplay.SelectOrientationVectors = 'None'
    sphDisplay.ScaleFactor = 0.2
    sphDisplay.SelectScaleArray = 'None'
    sphDisplay.GlyphType = 'Arrow'
    sphDisplay.GlyphTableIndexArray = 'None'
    sphDisplay.DataAxesGrid = 'GridAxesRepresentation'
    sphDisplay.PolarAxes = 'PolarAxesRepresentation'

    # change solid color
    sphDisplay.DiffuseColor = [1.0, 1.0, 0.5]







#------------------
# z axis
#------------------
cylinderZ = Cylinder()
cylinderZ.Radius = 0.05
cylinderZ.Center = [0.0, 15./2., 0.0]
cylinderZ.Height = h
cylinderZDisplay = Show(cylinderZ, renderView1)
cylinderZDisplay.ColorArrayName = [None, '']
cylinderZDisplay.DiffuseColor = [0.0, 1.0, 0.0]
# Default Cylinder is orientated along y-axis.
# To get z cylinder, rotate by 90 around x-axis.
cylinderZDisplay.Orientation = [90.0, 0.0, 0.0]
# cone z
coneZ = Cone()
# Properties modified on coneZ
coneZ.Resolution = 30
coneZ.Radius = 0.2
coneZ.Height = 0.4
coneZ.Center = [0.0, 0.0, h]
coneZ.Direction = [0., 0., 1.]
# show data in view
coneZDisplay = Show(coneZ, renderView1)
# trace defaults for the display properties.
coneZDisplay.Representation = 'Surface'
coneZDisplay.ColorArrayName = [None, '']
coneZDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coneZDisplay.SelectOrientationVectors = 'None'
coneZDisplay.ScaleFactor = 0.04000000357627869
coneZDisplay.SelectScaleArray = 'None'
coneZDisplay.GlyphType = 'Arrow'
coneZDisplay.GlyphTableIndexArray = 'None'
coneZDisplay.DataAxesGrid = 'GridAxesRepresentation'
coneZDisplay.PolarAxes = 'PolarAxesRepresentation'
# change solid color
coneZDisplay.DiffuseColor = [0.0, 1.0, 0.0]
for i in range(int(h)-1):
    # create a new 'Sphere'
    sph = Sphere()

    # Properties modified on sph
    sph.Center = [0., 0., i+1]
    sph.Radius = 0.2
    sph.ThetaResolution = 10
    sph.PhiResolution = 10

    # show data in view
    sphDisplay = Show(sph, renderView1)
    # trace defaults for the display properties.
    sphDisplay.Representation = 'Surface'
    sphDisplay.ColorArrayName = [None, '']
    sphDisplay.OSPRayScaleArray = 'Normals'
    sphDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    sphDisplay.SelectOrientationVectors = 'None'
    sphDisplay.ScaleFactor = 0.2
    sphDisplay.SelectScaleArray = 'None'
    sphDisplay.GlyphType = 'Arrow'
    sphDisplay.GlyphTableIndexArray = 'None'
    sphDisplay.DataAxesGrid = 'GridAxesRepresentation'
    sphDisplay.PolarAxes = 'PolarAxesRepresentation'

    # change solid color
    sphDisplay.DiffuseColor = [0.0, 1.0, 0.0]








deg=np.pi/180.
Mx,My,Mz = Mdipole
alpha = np.arccos(np.sqrt(Mx**2+My**2))
gamma = np.arctan(-Mx/My)
# cylinder M
cylinderM = Cylinder()
cylinderM.Radius = 0.05
cylinderM.Center = [0.0, h/2., 0.0]
cylinderM.Height = h
cylinderMDisplay = Show(cylinderM, renderView1)
cylinderMDisplay.ColorArrayName = [None, '']
cylinderMDisplay.DiffuseColor = [0.0, 0.0, 1.0]
# Default Cylinder is orientated along y-axis.
# To get z cylinder, rotate by 90 around x-axis.
cylinderMDisplay.Orientation = [alpha/deg, 1./137., gamma/deg]  #1./137. arbitrary
# cone z
coneM = Cone()
# Properties modified on coneM
coneM.Resolution = 30
coneM.Radius = 0.2
coneM.Height = 0.4
coneM.Center = h*Mdipole
coneM.Direction = Mdipole
# show data in view
coneMDisplay = Show(coneM, renderView1)
# trace defaults for the display properties.
coneMDisplay.Representation = 'Surface'
coneMDisplay.ColorArrayName = [None, '']
coneMDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coneMDisplay.SelectOrientationVectors = 'None'
coneMDisplay.ScaleFactor = 0.04000000357627869
coneMDisplay.SelectScaleArray = 'None'
coneMDisplay.GlyphType = 'Arrow'
coneMDisplay.GlyphTableIndexArray = 'None'
coneMDisplay.DataAxesGrid = 'GridAxesRepresentation'
coneMDisplay.PolarAxes = 'PolarAxesRepresentation'
# change solid color
coneMDisplay.DiffuseColor = [0.0, 0., 1.]

files = os.listdir(conf["run_path_derived"] + subdir)
for f in files:
    # create a new 'Legacy VTK Reader'
    if 'line' in f:
        field_linevtk = LegacyVTKReader(FileNames=[conf["run_path_derived"] + subdir + f])
        # show data in view
        field_linevtkDisplay = Show(field_linevtk, renderView1)
        # trace defaults for the display properties.
        field_linevtkDisplay.Representation = 'Surface'
        field_linevtkDisplay.ColorArrayName = [None, '']
        field_linevtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
        field_linevtkDisplay.SelectOrientationVectors = 'None'
        field_linevtkDisplay.ScaleFactor = 0.20896326303482057
        field_linevtkDisplay.SelectScaleArray = 'None'
        field_linevtkDisplay.GlyphType = 'Arrow'
        field_linevtkDisplay.GlyphTableIndexArray = 'None'
        field_linevtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
        field_linevtkDisplay.PolarAxes = 'PolarAxesRepresentation'


        # create a new 'Tube'
        tube1 = Tube(Input=field_linevtk, guiName=f)
        tube1.Scalars = [None, '']
        tube1.Vectors = [None, '1']
        tube1.Radius = 0.05

        # Properties modified on tube1
        tube1.Vectors = [None, '']

        # show data in view
        tube1Display = Show(tube1, renderView1)
        # trace defaults for the display properties.
        tube1Display.Representation = 'Surface'
        tube1Display.ColorArrayName = [None, '']
        tube1Display.OSPRayScaleArray = 'TubeNormals'
        tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        tube1Display.SelectOrientationVectors = 'None'
        tube1Display.ScaleFactor = 0.2129082262516022
        tube1Display.SelectScaleArray = 'None'
        tube1Display.GlyphType = 'Arrow'
        tube1Display.GlyphTableIndexArray = 'None'
        tube1Display.DataAxesGrid = 'GridAxesRepresentation'
        tube1Display.PolarAxes = 'PolarAxesRepresentation'

        # hide data in view
        Hide(field_linevtk, renderView1)

        # change solid color
        if 'event' in f:
            tube1Display.DiffuseColor = [0.8862745098039215, 0.0, 0.0]
        elif 'longitude' in f:
            tube1Display.DiffuseColor = [0., 0.0, 1.]
        elif 'B_' in f:
            tube1Display.DiffuseColor = [0., 1.0, 0.]
        elif 'J_' in f:
            tube1Display.DiffuseColor = [0., 0., 1.]
        else:
            tube1Display.DiffuseColor = [0.0, 0., 0.0]

renderView1.Update()

# find source
slice1 = FindSource('Slice1')

# rename source object
RenameSource('xz_plane', slice1)

# find source
slice2 = FindSource('Slice2')

# rename source object
RenameSource('xy_plane', slice2)

# find source
slice3 = FindSource('Slice3')

# rename source object
RenameSource('yz_plane', slice3)

# find source
slice4 = FindSource('Slice4')

# rename source object
RenameSource('cut_plane', slice4)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 0
renderView1.OrientationAxesVisibility = 0

point_scalarsLUT = GetColorTransferFunction(var)
point_scalarsLUT.RescaleTransferFunction(-0.0005, 0.0005)

renderView1.CameraPosition = [7.1183569, 56.134761, 6.60090036]
renderView1.CameraFocalPoint = [-0.445446427, 0.115471756, -0.440184]
renderView1.CameraViewUp = [-0.05222445, -0.117594, 0.9916875635]
renderView1.CenterOfRotation = [0., 0., 0.]
renderView1.RotationFactor = 1
renderView1.CameraViewAngle = 30
renderView1.CameraParallelScale = 108.426242211007
renderView1.CameraParallelProjection = 0

'''
<?xml version="1.0"?>
<PVCameraConfiguration description="ParaView camera configuration" version="1.0">
  <Proxy group="views" type="RenderView" id="2076" servers="21">
    <Property name="CameraPosition" id="2076.CameraPosition" number_of_elements="3">
      <Element index="0" value="7.11835691704117"/>
      <Element index="1" value="56.1347613281343"/>
      <Element index="2" value="6.60090036368214"/>
    </Property>
    <Property name="CameraFocalPoint" id="2076.CameraFocalPoint" number_of_elements="3">
      <Element index="0" value="-0.445446427868706"/>
      <Element index="1" value="0.115471756220899"/>
      <Element index="2" value="-0.440184037664585"/>
    </Property>
    <Property name="CameraViewUp" id="2076.CameraViewUp" number_of_elements="3">
      <Element index="0" value="-0.0522244587063913"/>
      <Element index="1" value="-0.117594142115094"/>
      <Element index="2" value="0.991687563526457"/>
    </Property>
    <Property name="CenterOfRotation" id="2076.CenterOfRotation" number_of_elements="3">
      <Element index="0" value="-0.445446427868706"/>
      <Element index="1" value="0.115471756220899"/>
      <Element index="2" value="-0.440184037664585"/>
    </Property>
    <Property name="RotationFactor" id="2076.RotationFactor" number_of_elements="1">
      <Element index="0" value="1"/>
    </Property>
    <Property name="CameraViewAngle" id="2076.CameraViewAngle" number_of_elements="1">
      <Element index="0" value="30"/>
    </Property>
    <Property name="CameraParallelScale" id="2076.CameraParallelScale" number_of_elements="1">
      <Element index="0" value="108.426242211007"/>
    </Property>
    <Property name="CameraParallelProjection" id="2076.CameraParallelProjection" number_of_elements="1">
      <Element index="0" value="0"/>
      <Domain name="bool" id="2076.CameraParallelProjection.bool"/>
    </Property>
  </Proxy>
</PVCameraConfiguration>
'''
