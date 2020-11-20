# adapted from https://gitlab.kitware.com/vtk/vtk/blob/32ae6c76b70c8da05cec8cabe3ccda8988388561/Filters/General/Testing/Python/streamTracer.py
# see ~/home/gary/VTKData/streamTracer.py

import vtk
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa


def trace(X, view=False, write=False):
    datafile = "/tmp/tmp.vtk"

    # read data
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(datafile)
    reader.Update()
    #force a read to occur


    rk = vtk.vtkRungeKutta45()
    # Create source for streamtubes
    streamer = vtk.vtkStreamTracer()
    streamer.SetInputConnection(reader.GetOutputPort())
    streamer.SetStartPosition(X[0],X[1],X[2])#SetSourceData()
    streamer.SetMaximumPropagation(20) ###
    #streamer.SetIntegrationStepUnit(2) # apperars overiden by next lines, see https://vtk.org/doc/nightly/html/classvtkStreamTracer.html#afe365e81e110f354065f5adc8401d589
    streamer.SetMinimumIntegrationStep(0.01)
    streamer.SetMaximumIntegrationStep(0.1)
    streamer.SetInitialIntegrationStep(0.01)
    #streamer.SetIntegrationDirection()
    #streamer.SetIntegrationDirectionToForward()
    #streamer.SetIntegrationDirectionToBackward()
    streamer.SetIntegrationDirectionToBoth()
    streamer.SetIntegrator(rk)
    streamer.SetRotationScale(0.5)
    streamer.SetMaximumError(1.0e-8)

    print('#####################')
    print(streamer)
    print(type(streamer))
    print('#####################')
    print(streamer.GetOutputPort())
    print('#####################')

    if write:
        writer = vtk.vtkDataSetWriter()
        writer.SetInputConnection(streamer.GetOutputPort())
        writer.SetFileName("out.vtk")
        writer.Write()

    if view:
        outline = vtk.vtkStructuredGridOutlineFilter()
        outline.SetInputConnection(reader.GetOutputPort())
        mapOutline = vtk.vtkPolyDataMapper()
        mapOutline.SetInputConnection(outline.GetOutputPort())
        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(mapOutline)
        outlineActor.GetProperty().SetColor(0,0,0)

        mapStream = vtk.vtkPolyDataMapper()
        print(streamer.GetOutputPort())
        print(type(streamer.GetOutputPort()))
        mapStream.SetInputConnection(streamer.GetOutputPort())
        mapStream.SetScalarRange(reader.GetOutput().GetScalarRange())
        streamActor = vtk.vtkActor()
        streamActor.SetMapper(mapStream)

        return [outlineActor, streamActor, streamer]


    #https://stackoverflow.com/questions/38504907/reading-a-vtk-polydata-file-and-converting-it-into-numpy-array

    ## either order works ##
    polydata = streamer.GetOutput()
    streamer.Update()

    return dsa.WrapDataObject(polydata).Points



    #https://searchcode.com/file/25190165/Filters/Core/Testing/Python/fieldToRGrid.py/
    #ds2do = vtk.vtkDataSetToDataObjectFilter()
    #ds2do.SetInputConnection(streamer.GetOutputPort())
    #writer = vtk.vtkDataObjectWriter()
    #writer.SetInputConnection(ds2do.GetOutputPort())

    #ds_att = vtk.vtkDataSetAttributes()
    #ds_att.SetInputConnection(streamer.GetOutputPort())

    #pnts = vtk.vtkPoints()
    #pnts.SetData(streamer.GetOutputPort())
    #print(pnts)

    #rdr = vtk.vtkDataSetReader()
    #rdr.SetInputConnection(streamer.GetOutputPort())

def main():
    ren1 = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    # this test has some wireframe geometry
    # Make sure multisampling is disabled to avoid generating multiple
    # regression images
    # renWin SetMultiSamples 0
    renWin.AddRenderer(ren1)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    X = np.array([0.1,2.1,0.5])
    outlineActor, streamActor, streamer = trace(X, view=True)

    #print(streamer.GetOutputPort())

    ren1.AddActor(outlineActor)
    ren1.AddActor(streamActor)

    ren1.SetBackground(0.4,0.4,0.5)
    cam = ren1.GetActiveCamera()
    cam.SetPosition(20.,20.,20.) ###
    cam.SetFocalPoint(0.,0.,0.) ###
    cam.SetViewUp(0.311311,0.279912,0.908149)
    cam.SetClippingRange(1.12294,16.6226)
    renWin.SetSize(900,600)
    iren.Initialize()

    ### added ###
    renWin.Render()
    iren.Start()

    #print(streamer.GetOutputPort())

#vtkAlgorithmOutput

if __name__=='__main__':
    main()

    X = np.array([0.1,2.1,0.5])
    arr = trace(X, view=False, write=True)
    print(arr)
