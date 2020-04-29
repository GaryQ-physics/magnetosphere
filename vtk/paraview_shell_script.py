# paraview_shell_script

#    execfile('paraview_shell_script.py')   #copy to shell


#fieldFile_path = m_path + 'magnetosphere/vtk/kameleon_field_paraview.vtk'

import sys
import os
import numpy as np
#print os.path.dirname(os.path.abspath(__file__)) + '/../'
#sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../' )
sys.path.append( './../' )  # for some reason the other way wont work in paraview
from config_paths import config
conf = config()
sys.path.append(conf["m_path"] + 'magnetosphere/events/')
import pos_sun as ps
from cut_plane import U1,Mdipole,U2,U3 #note cut plane has to be run seperately once first
#from paraview.simple import *              already in shell automatically
Nlong=5
Nb = 6
N=Nb+1+Nlong

# create a new 'Legacy VTK Reader'
kameleon_structured_gridvtk = LegacyVTKReader(FileNames=[conf["m_path"] + 'magnetosphere/vtk/kameleon_structured_grid.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [932, 802]

# show data in view
kameleon_structured_gridvtkDisplay = Show(kameleon_structured_gridvtk, renderView1)
# trace defaults for the display properties.
kameleon_structured_gridvtkDisplay.Representation = 'Outline'
kameleon_structured_gridvtkDisplay.ColorArrayName = ['POINTS', '']
kameleon_structured_gridvtkDisplay.OSPRayScaleArray = 'point_scalars'
kameleon_structured_gridvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
kameleon_structured_gridvtkDisplay.SelectOrientationVectors = 'None'
kameleon_structured_gridvtkDisplay.ScaleFactor = 21.0
kameleon_structured_gridvtkDisplay.SelectScaleArray = 'point_scalars'
kameleon_structured_gridvtkDisplay.GlyphType = 'Arrow'
kameleon_structured_gridvtkDisplay.GlyphTableIndexArray = 'point_scalars'
kameleon_structured_gridvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
kameleon_structured_gridvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
kameleon_structured_gridvtkDisplay.ScalarOpacityUnitDistance = 5.766431907004487

# reset view to fit data
renderView1.ResetCamera()

# change representation type
kameleon_structured_gridvtkDisplay.SetRepresentationType('Points')

# create a new 'Slice'
slice1 = Slice(Input=kameleon_structured_gridvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-95.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# get color transfer function/color map for 'point_scalars'
point_scalarsLUT = GetColorTransferFunction('point_scalars')

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'point_scalars']
slice1Display.LookupTable = point_scalarsLUT
slice1Display.OSPRayScaleArray = 'point_scalars'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 2.0
slice1Display.SelectScaleArray = 'point_scalars'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'point_scalars'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)


# create a new 'Slice'
slice2 = Slice(Input=kameleon_structured_gridvtk)
slice2.SliceType = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.0, 0.0, 0.0]
slice2.SliceType.Normal = U3

# get color transfer function/color map for 'point_scalars'
point_scalarsLUT = GetColorTransferFunction('point_scalars')

# show data in view
slice2Display = Show(slice2, renderView1)
# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['POINTS', 'point_scalars']
slice2Display.LookupTable = point_scalarsLUT
slice2Display.OSPRayScaleArray = 'point_scalars'
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 2.0
slice2Display.SelectScaleArray = 'point_scalars'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'point_scalars'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# show color bar/color legend
slice2Display.SetScalarBarVisibility(renderView1, True)


# hide data in view
Hide(kameleon_structured_gridvtk, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# create a new 'Sphere'
sphere1 = Sphere()

# Properties modified on sphere1
sphere1.Center = [0.01, 0.0, 0.0]
sphere1.Radius = 1.0
sphere1.ThetaResolution = 20
sphere1.PhiResolution = 20

# show data in view
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.Representation = 'Surface'
sphere1Display.ColorArrayName = [None, '']
sphere1Display.OSPRayScaleArray = 'Normals'
sphere1Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere1Display.SelectOrientationVectors = 'None'
sphere1Display.ScaleFactor = 0.2
sphere1Display.SelectScaleArray = 'None'
sphere1Display.GlyphType = 'Arrow'
sphere1Display.GlyphTableIndexArray = 'None'
sphere1Display.DataAxesGrid = 'GridAxesRepresentation'
sphere1Display.PolarAxes = 'PolarAxesRepresentation'

# find source
#glyph1 = FindSource('Glyph1')

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Sphere'
sphere2 = Sphere()

# Properties modified on sphere2
sphere2.Center = [-0.01, 0.0, 0.0]
sphere2.Radius = 1.0
sphere2.ThetaResolution = 20
sphere2.PhiResolution = 20

# show data in view
sphere2Display = Show(sphere2, renderView1)
# trace defaults for the display properties.
sphere2Display.Representation = 'Surface'
sphere2Display.ColorArrayName = [None, '']
sphere2Display.OSPRayScaleArray = 'Normals'
sphere2Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere2Display.SelectOrientationVectors = 'None'
sphere2Display.ScaleFactor = 0.2
sphere2Display.SelectScaleArray = 'None'
sphere2Display.GlyphType = 'Arrow'
sphere2Display.GlyphTableIndexArray = 'None'
sphere2Display.DataAxesGrid = 'GridAxesRepresentation'
sphere2Display.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# change solid color
sphere2Display.DiffuseColor = [0.0, 0.0, 0.0]

# cylinder x
cylinderX = Cylinder()
cylinderX.Radius = 0.02
cylinderX.Center = [0.0, 1.0, 0.0]
cylinderX.Height = 2.0
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
coneX.Center = [2.0, 0.0, 0.0]
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

# y
cylinderY = Cylinder()
cylinderY.Radius = 0.02
cylinderY.Center = [0.0, 1., 0.0]
cylinderY.Height = 2.0
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
coneY.Center = [0.0, 2.0, 0.0]
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


# z
cylinderZ = Cylinder()
cylinderZ.Radius = 0.02
cylinderZ.Center = [0.0, 1., 0.0]
cylinderZ.Height = 2.0
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
coneZ.Center = [0.0, 0.0, 2.0]
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


deg=np.pi/180.
Mx,My,Mz = Mdipole
alpha = np.arccos(np.sqrt(Mx**2+My**2))
gamma = np.arctan(-Mx/My)
# cylinder M
cylinderM = Cylinder()
cylinderM.Radius = 0.02
cylinderM.Center = [0.0, 1., 0.0]
cylinderM.Height = 2.0
cylinderMDisplay = Show(cylinderM, renderView1)
cylinderMDisplay.ColorArrayName = [None, '']
cylinderMDisplay.DiffuseColor = [0.0, 0.0, 1.0]
# Default Cylinder is orientated along y-axis.
# To get z cylinder, rotate by 90 around x-axis.
cylinderMDisplay.Orientation = [alpha/deg, 3.14, gamma/deg]  #3.14 arbitrary
# cone z
coneM = Cone()
# Properties modified on coneM
coneM.Resolution = 30
coneM.Radius = 0.2
coneM.Height = 0.4
coneM.Center = 2*Mdipole
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

for i in range(N):
    # create a new 'Legacy VTK Reader'
    field_linevtk = LegacyVTKReader(FileNames=[conf["m_path"] + 'magnetosphere/vtk/' + 'field_line'+str(i)+'.vtk'])

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
    tube1 = Tube(Input=field_linevtk)
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
    if(i==0):
        tube1Display.DiffuseColor = [0.8862745098039215, 0.0, 0.0]
    elif(i > Nb):
        tube1Display.DiffuseColor = [0., 0.0, 1.]
    else:
        tube1Display.DiffuseColor = [0.0, 0.88, 0.0]

renderView1.Update()

