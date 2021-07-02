# VTK file types

There are two main VTK file types: Legacy and XML.

see:
- `https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf`
- `https://vtk.org/Wiki/VTK_XML_Formats`

## Legacy

Filetype extension: `.vtk`

## XML

Filetype extentions: `.vtu` (for unstructured grids), `.vts` (for structured grids), and others.
Note there are different extenstions for different datasets 


# Libraries

## vtk (python2 and python3)

A python library that minimally wraps the VTK c++ library.
As such, it utilizes the native VTK "pipeline style" of doing things, and has the full functionality.
It can read (and write?) files, as well as perform all the interpolation, streamlines, ect that VTK can do.

## meshio (python3)

### info
A python library to convert between many different file types that can decribe meshes,
including Legacy and XML VTK files.

### notes
Since its not vtk specific, it imports vtk structured grids as unstructured grids.
So if you import a structured grid vtk file and export it again, it will be an
unstructured grid vtk file. Legacy vtk files exported are version 4.2

## pyvista (python3)

### info
A high level library for alot of general 3D visualization and mesh analysis stuff.
It has meshio and vtk as dependency, among other things.

## pyvtk (python2 and python3)

### info
A python library specifically for reading and writing VTK files.

### notes
Would be great, but its almost prohibitively slow.

## uvw

Supposedly a small utility library to write XML VTK files from data contained in Numpy arrays.
I never got it to work however.

## vtk\_export \[part of magnetovis\] (python2 and python3)

### info
A single python function, independent of any packages but numpy,
that writes Legacy VTK files based on python data and numpy arrays passed to it.
There are no custom mesh objects and no reading files capability.

### notes
This is the one we wrote ourselves.
It only has only the features we needed at the time.
It is very fast.
