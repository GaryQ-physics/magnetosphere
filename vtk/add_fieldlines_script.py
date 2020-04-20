# add_fieldlines_script

Nb = 6
N=Nb+1

for i in range(N):
    # create a new 'Legacy VTK Reader'
    field_linevtk = LegacyVTKReader(FileNames=['field_line'+str(i)+'.vtk'])

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
    else:
        tube1Display.DiffuseColor = [0.0, 0.88, 0.0]

renderView1.Update()
