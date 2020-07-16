"""
Usage:
    Start Paraview from this directroy and in Paraview Python shell, enter
        execfile('axes_paraview.py')
    or on command line enter
        paraview5.8 --script=axes_paraview.py
"""

from paraview.simple import *  # needed for paraview 5.8+

renderView1 = GetActiveViewOrCreate('RenderView')

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

# find source
cylinder1 = FindSource('Cylinder1')
cylinder2 = FindSource('Cylinder2')
cylinder3 = FindSource('Cylinder3')

# rename source object
RenameSource('x-axis_rod', cylinder1)
RenameSource('y-axis_rod', cylinder2)
RenameSource('z-axis_rod', cylinder3)

