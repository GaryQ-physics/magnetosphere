import os
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

def printl(obj):
    print('\n\n')
    print(type(obj))
    print(obj)
    print('\n\n')

fname = '/home/gary/temp/'+'3d__var_3_e20031120-070000-000.vtu'

fname = fname[:-4]+'.vtk' # converted from vtu to vtk using meshio python library 
                          # careful using paraview to save loaded vtu, as it can remove the cells

#fname = '/tmp/tmp.vtk'
print(fname+'\n\n')

#################################################
UseNewPipeline = False # https://vtk.org/Wiki/VTK/Tutorials/New_Pipeline
debug = False


## open file
if not os.path.exists(fname): raise FileNotFoundError ('no file ' + fname)
extention = fname[-4:]
if extention == '.vtk':
    reader = vtk.vtkGenericDataObjectReader()
elif extention == '.vtu':
    reader = vtk.vtkXMLUnstructuredGridReader() # https://stackoverflow.com/questions/54044958/reading-data-from-a-raw-vtk-vtu-file
reader.SetFileName(fname)
reader.Update() # forces loading the file into memory
if debug: printl(reader)
# get nessesary vtk info from Reader object
if UseNewPipeline:
    output_port = reader.GetOutputPort()
    printl(output_port)
    # how to set which field to use?
else:
    output = reader.GetOutput()
    if debug: printl(output)
    if True:
        output.GetPointData().SetActiveAttribute('b', 1)
        # https://vtk.org/doc/nightly/html/classvtkDataSetAttributes.html#ac0d39efa092df157d210578a5982f7e8
        # Attribute types are: vtkDataSetAttributes::SCALARS = 0 vtkDataSetAttributes::VECTORS = 1 vtkDataSetAttributes::NORMALS = 2 vtkDataSetAttributes::TCOORDS = 3 vtkDataSetAttributes::TENSORS = 4 vtkDataSetAttributes::GLOBALIDS = 5 vtkDataSetAttributes::PEDIGREEIDS = 6 vtkDataSetAttributes::EDGEFLAG = 7 vtkDataSetAttributes::TANGENTS = 8 Returns the index of the array if successful, -1 if the array is not in the list of arrays. 
    else:
        output.GetPointData().SetActiveVectors('b')
    if debug: printl(output)


## feed the above vtk info into the nesseary streamline integration objects
# Create integrator
rk = vtk.vtkRungeKutta45()
# Create source for streamtubes
streamer = vtk.vtkStreamTracer()
if UseNewPipeline:
    streamer.SetInputConnection(output_port)
else:
    streamer.SetInputDataObject(output)
IC = np.array([1.,1.,1.])
streamer.SetStartPosition(IC)
#streamer.SetStartPosition(1., 1., 1.)# also works, but either way cannot pass multiple IC's in an array
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

#https://stackoverflow.com/questions/38504907/reading-a-vtk-polydata-file-and-converting-it-into-numpy-array
## either order works ##
polydata = streamer.GetOutput()

if debug: printl(polydata)
if debug: printl(streamer)

streamer.Update() # forces the computation of the stream lines

if debug: printl(polydata)
if debug: printl(streamer)

arr = dsa.WrapDataObject(polydata).Points  # convert the result to array

print(arr.shape)
print(arr)
