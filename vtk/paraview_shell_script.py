# paraview_shell_script

#    execfile('paraview_shell_script.py')   #copy to shell

k_path = '/home/gary/magnetosphere/'
fieldFile_path = k_path + 'vtk/kameleon_field_paraview.vtk'

import sys
import numpy as np
sys.path.append(k_path + 'events/')
import pos_sun as ps
from cut_plane import U1,Mdipole,U2,U3


#### import the simple module from the paraview
#from paraview.simple import *              already in shell automatically

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
kameleon_field_paraviewvtk = LegacyVTKReader(FileNames=[fieldFile_path])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [790, 802]

# show data in view
kameleon_field_paraviewvtkDisplay = Show(kameleon_field_paraviewvtk, renderView1)
# trace defaults for the display properties.
kameleon_field_paraviewvtkDisplay.Representation = 'Surface'
kameleon_field_paraviewvtkDisplay.ColorArrayName = [None, '']
kameleon_field_paraviewvtkDisplay.OSPRayScaleArray = 'point_vectors'
kameleon_field_paraviewvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
kameleon_field_paraviewvtkDisplay.SelectOrientationVectors = 'point_vectors'
kameleon_field_paraviewvtkDisplay.SelectScaleArray = 'None'
kameleon_field_paraviewvtkDisplay.GlyphType = 'Arrow'
kameleon_field_paraviewvtkDisplay.GlyphTableIndexArray = 'None'
kameleon_field_paraviewvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
kameleon_field_paraviewvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
kameleon_field_paraviewvtkDisplay.SetRepresentationType('Points')

# create a new 'Glyph'
glyph1 = Glyph(Input=kameleon_field_paraviewvtk,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'point_vectors']
glyph1.GlyphTransform = 'Transform2'

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'GlyphVector'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 1.1981441497802734
glyph1Display.SelectScaleArray = 'GlyphVector'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'GlyphVector'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# PLANE 1
plane1 = Plane()
w0=-4*U1-4*U2
w1=-4*U1+4*U2
w2=50*U1-4*U2
# Properties modified on plane1
plane1.Origin = w0
plane1.Point1 = w1
plane1.Point2 = w2
#plane1.Origin = [-2.0, -1.0, 0.0]
#plane1.Point1 = [1.0, 2.0, 3.0]
#plane1.Point2 = [4.0, -3.5, 1.0]
# show data in view
plane1Display = Show(plane1, renderView1)
# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.ColorArrayName = [None, '']
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 0.9500000000000001
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'
# change solid color
plane1Display.DiffuseColor = [0.6666666666666666, 0.0, 0.4980392156862745]
# Properties modified on plane1Display
plane1Display.Opacity = 0.5

# PLANE 2
plane2 = Plane()
w0=[4., 0., 4.]
w1=[4., 0., -4.]
w2=[-200., 0., 4.]
# Properties modified on plane1
plane2.Origin = w0
plane2.Point1 = w1
plane2.Point2 = w2
# show data in view
plane2Display = Show(plane2, renderView1)
# trace defaults for the display properties.
plane2Display.Representation = 'Surface'
plane2Display.ColorArrayName = [None, '']
plane2Display.OSPRayScaleArray = 'Normals'
plane2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane2Display.SelectOrientationVectors = 'None'
plane2Display.ScaleFactor = 0.9500000000000001
plane2Display.SelectScaleArray = 'None'
plane2Display.GlyphType = 'Arrow'
plane2Display.GlyphTableIndexArray = 'None'
plane2Display.DataAxesGrid = 'GridAxesRepresentation'
plane2Display.PolarAxes = 'PolarAxesRepresentation'
# change solid color
plane2Display.DiffuseColor = [1.0, 0.6666666666666666, 0.0]
# Properties modified on plane1Display
plane2Display.Opacity = 0.5

# PLANE 3
plane3 = Plane()
w0=[4., 4., 0.]
w1=[4., -4., 0.]
w2=[-200., 4., 0.]
# Properties modified on plane1
plane3.Origin = w0
plane3.Point1 = w1
plane3.Point2 = w2
# show data in view
plane3Display = Show(plane3, renderView1)
# trace defaults for the display properties.
plane3Display.Representation = 'Surface'
plane3Display.ColorArrayName = [None, '']
plane3Display.OSPRayScaleArray = 'Normals'
plane3Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane3Display.SelectOrientationVectors = 'None'
plane3Display.ScaleFactor = 0.9500000000000001
plane3Display.SelectScaleArray = 'None'
plane3Display.GlyphType = 'Arrow'
plane3Display.GlyphTableIndexArray = 'None'
plane3Display.DataAxesGrid = 'GridAxesRepresentation'
plane3Display.PolarAxes = 'PolarAxesRepresentation'
# change solid color
plane3Display.DiffuseColor = [1.0, 0.6666666666666666, 0.0]
# Properties modified on plane1Display
plane3Display.Opacity = 0.5

# PLANE 4
plane2 = Plane()
w0=[0., 4., 4.]
w1=[0., 4., -4.]
w2=[0., -4., 4.]
# Properties modified on plane1
plane2.Origin = w0
plane2.Point1 = w1
plane2.Point2 = w2
# show data in view
plane2Display = Show(plane2, renderView1)
# trace defaults for the display properties.
plane2Display.Representation = 'Surface'
plane2Display.ColorArrayName = [None, '']
plane2Display.OSPRayScaleArray = 'Normals'
plane2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane2Display.SelectOrientationVectors = 'None'
plane2Display.ScaleFactor = 0.9500000000000001
plane2Display.SelectScaleArray = 'None'
plane2Display.GlyphType = 'Arrow'
plane2Display.GlyphTableIndexArray = 'None'
plane2Display.DataAxesGrid = 'GridAxesRepresentation'
plane2Display.PolarAxes = 'PolarAxesRepresentation'
# change solid color
plane2Display.DiffuseColor = [1.0, 0.6666666666666666, 0.0]
# Properties modified on plane1Display
plane2Display.Opacity = 0.5


# update the view to ensure updated data information
renderView1.Update()

# get active view
#renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 802]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# find source
legacyVTKReader1 = FindSource('LegacyVTKReader1')

# find source
glyph1 = FindSource('Glyph1')

# change color vector field arrows
glyph1Display.DiffuseColor = [0.0, 0.0, 0.8705882352941177]

# create a new 'Sphere'
sphere1 = Sphere()

# find source
plane1 = FindSource('Plane1')

# Properties modified on plane1
plane1.Origin = [4.32941593011924, -2.80002404457948, -2.32723506629789]
plane1.Point1 = [0.713399586162049, -2.89333858275426, 4.80829001579661]
plane1.Point2 = [-29.7095888047794, 35.6301736899233, -19.0743559754143]

# find source
legacyVTKReader1 = FindSource('LegacyVTKReader1')

# Properties modified on sphere1
sphere1.Center = [0.01, 0.0, 0.0]
sphere1.Radius = 1.0
sphere1.ThetaResolution = 20
sphere1.PhiResolution = 20

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 802]

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
plane2 = FindSource('Plane2')

# find source
glyph1 = FindSource('Glyph1')

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


# update the view to ensure updated data information
renderView1.Update()
