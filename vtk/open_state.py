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
LoadState('/home/gary/magnetosphere/vtk/test_state.pvsm', LoadStateDataFileOptions='Use File Names From State',
    DataDirectory='/home/gary/magnetosphere/vtk',
    OnlyUseFilesInDataDirectory=0,
    LegacyVTKReader1FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/kameleon_structured_grid_dBy_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader2FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line0_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader3FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line1_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader4FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line2_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader5FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line3_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader6FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line4_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader7FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line5_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader8FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line6_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader9FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line7_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader10FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line8_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader11FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line9_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader12FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line10_2003:11:20T07:00:00.vtk'],
    LegacyVTKReader13FileNames=['/home/gary/magnetosphere/data/SCARR5_GM_IO2-derived/field_line11_2003:11:20T07:00:00.vtk'])
