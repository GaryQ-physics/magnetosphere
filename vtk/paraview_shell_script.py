# paraview_shell_script

#    execfile('paraview_shell_script.py')   #copy to shell

k_path = '/home/gary/magnetosphere/'
fieldFile_path = k_path + 'vtk/kameleon_field_paraview.vtk'

import sys
import numpy as np
sys.path.append(k_path + 'events/')
from cut_plane import U1,U2,U3


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


# create a new 'Plane'
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

# find source
plane1 = FindSource('Plane1')


# change solid color
glyph1Display.DiffuseColor = [0.0, 0.0, 0.8705882352941177]



# get active source.
#plane1 = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 800]

# get display properties
plane1Display = GetDisplayProperties(plane1, view=renderView1)

# change solid color
plane1Display.DiffuseColor = [0.6666666666666666, 0.0, 0.4980392156862745]

# Properties modified on plane1Display
plane1Display.Opacity = 0.5

# get active source.
#plane2 = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [964, 800]

# get display properties
plane2Display = GetDisplayProperties(plane2, view=renderView1)

# change solid color
plane2Display.DiffuseColor = [1.0, 0.6666666666666666, 0.0]

# Properties modified on plane1Display
plane2Display.Opacity = 0.5

# change solid color
plane1Display.DiffuseColor = [0.6666666666666666, 0.0, 0.4980392156862745]

# Properties modified on plane1Display
plane1Display.Opacity = 0.5

# update the view to ensure updated data information
renderView1.Update()
