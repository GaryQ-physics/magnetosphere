import paraview.simple as pvs

# https://trac.version.fz-juelich.de/vis/wiki/Examples/Ear5Animating

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

textureMaptoSphere = TextureMaptoSphere(Input=sphere1)
textureMaptoSphere.PreventSeam = 0
textureMaptoSphereDisplay = Show(textureMaptoSphere, renderView1)
textureMaptoSphereDisplay.Representation = 'Surface'
texProxy = servermanager.CreateProxy("textures", "ImageTexture")
texProxy.GetProperty("FileName").SetElement(0, "/Users/robertweigel/git/students/gquaresi/magnetosphere/data/topography/world.topo.200401.3x5400x2700.png-ParaView.png")
texProxy.UpdateVTKObjects()
textureMaptoSphereDisplay.Texture = texProxy

Render()