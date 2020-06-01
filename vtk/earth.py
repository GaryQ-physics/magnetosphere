import paraview.simple as pvs

# Based on
# https://trac.version.fz-juelich.de/vis/wiki/Examples/Ear5Animating

sys.path.append( './../' )
from config import conf
#from config_paths import config
#conf = config()

fname = '/home/gary/magnetosphere/vtk/rotated_sphere.vtk'

file = conf["data_path"] + "topography/world.topo.200401.3x5400x2700.png-ParaView.png"

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')


sphere1 = pvs.Sphere(guiName="Earth")

sphere1.Radius = 1.0
sphere1.ThetaResolution = 100
sphere1.PhiResolution = 100
sphere1.StartTheta = 1e-05

sphere1Display = Show(sphere1, renderView1)
sphere1Display.Representation = 'Surface'
sphere1Display.Opacity = 1.0


'''


# create a new 'Transform'
transform1 = Transform(Input=sphere1)
transform1.Transform = 'Transform'
# Properties modified on transform1.Transform
transform1.Transform.Rotate = [10.0, 20.0, -30.0]
# show data in view
transform1Display = Show(transform1, renderView1)
# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = [None, '']
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 0.199983811378479
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
# update the view to ensure updated data information
renderView1.Update()

'''

'''

# create a new 'Legacy VTK Reader'
rotated_spherevtk = LegacyVTKReader(FileNames=['/home/gary/magnetosphere/vtk/rotated_sphere.vtk'])
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

'''


textureMaptoSphere = TextureMaptoSphere(Input=sphere1, Point=[0,0,1])   #(Input=rotated_spherevtk)
textureMaptoSphere.PreventSeam = 0
textureMaptoSphereDisplay = Show(textureMaptoSphere, renderView1)
textureMaptoSphereDisplay.Representation = 'Surface'
texProxy = servermanager.CreateProxy("textures", "ImageTexture")
texProxy.GetProperty("FileName").SetElement(0, file)
texProxy.UpdateVTKObjects()
textureMaptoSphereDisplay.Texture = texProxy

Render()
