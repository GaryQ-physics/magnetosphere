# paraview_shell_script

#    execfile('paraview_shell_script.py')   #copy to shell

k_path = '/home/gary/magnetosphere/'
fieldFile_path = k_path + 'vtk/kameleon_field_paraview.vtk'
lineFile_path = k_path + 'vtk/field_line.vtk'

import sys
import numpy as np
sys.path.append(k_path + 'events/')
from cut_plane import U1
from cut_plane import U2
from cut_plane import U3


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

# create a new 'Plane'
plane1 = Plane()

w0=-4*U1-4*U2
w1=4*U1-4*U2
w2=4*U2-4*U1

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

# update the view to ensure updated data information
renderView1.Update()


# get active view
#renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 802]

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# create a new 'Legacy VTK Reader'
field_linevtk = LegacyVTKReader(FileNames=[lineFile_path])

# find source
legacyVTKReader1 = FindSource('LegacyVTKReader1')

# find source
glyph1 = FindSource('Glyph1')

# find source
plane1 = FindSource('Plane1')

# Properties modified on plane1
plane1.Origin = [4.32941593011924, -2.80002404457948, -2.32723506629789]
plane1.Point1 = [-0.713399586162049, 2.89333858275426, -4.80829001579661]
plane1.Point2 = [0.713399586162049, -2.89333858275426, 4.80829001579661]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 802]

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

# update the view to ensure updated data information
renderView1.Update()

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
tube1Display.DiffuseColor = [0.8862745098039215, 0.0, 0.0]

# change solid color
glyph1Display.DiffuseColor = [0.0, 0.0, 0.8705882352941177]

# update the view to ensure updated data information
renderView1.Update()


