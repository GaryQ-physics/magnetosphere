###copy to shell
'''
execfile('open_state.py')
'''
###

#### import the simple module from the paraview
#from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [989, 802]

# destroy renderView1
Delete(renderView1)
del renderView1

# load state
LoadState('../data/SCARR5_GM_IO2-derived/test_state.pvsm', 
    LoadStateDataFileOptions='Search files under specified directory',
    DataDirectory='../data/SCARR5_GM_IO2-derived/',
    OnlyUseFilesInDataDirectory=1,
    LegacyVTKReader1FileNames=['kameleon_structured_grid_dBy_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader2FileNames=['field_line0_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader3FileNames=['field_line1_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader4FileNames=['field_line2_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader5FileNames=['field_line3_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader6FileNames=['field_line4_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader7FileNames=['field_line5_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader8FileNames=['field_line6_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader9FileNames=['field_line7_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader10FileNames=['field_line8_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader11FileNames=['field_line9_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader12FileNames=['field_line10_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader13FileNames=['field_line11_2003:11:20T07:00:00.vtk'])
