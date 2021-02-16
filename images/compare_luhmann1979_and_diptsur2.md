# Introduction

Here we will look at plots from from the DIPTSUR2 swmf simulation run, 
and the LUHMANN1979 simple analytic magnetosphere model.

The DIPTSUR2 run
was run on a block tree grid, with points and values output to a
datafile (eventually converted to vtu).

The LUHMANN1979 is an analytic model from ![a paper](https://doi.org/10.1029/JA084iA08p04405)
which includes a callable function for magnetic field. This function was
then evaluated at each point of the DIPTSUR2 run grid, in order to make
a mock vtu file.

All the plots shown in the next section were done using points sampled from this native block tree grid.

The three swmf variables used here are: "b", "b1", "j". Here "b" represents
the full magnetic field, "b1" represent the magnetic field minus the contribution
from the earths dipole field, and "j" represents the current density.

In the DIPTSUR2 run, "j" is not actually part of the simulation equations. 
Suposedly "j" is calculated as the curl of "b" (multiplied by mu_0) evaluated
on the native block tree grid, using a discrete derivative formula. If this is
this case, then I would expect taking the descrete curl of "b" from the datafile
to exacly reporduce "j".

In the LUHMANN1979 mock run, the magnetic field model includes a dipole term,
and a non dipole term $ bla $. Here "b" was taken to be the full magnetic field and
"b1" the non dipole term, evaluated at the native grid points. Then, the analytic
magnetic field function was analytically differentiated to get the curl.
This analytic curl, multiplied by mu_0, was evalueated at the grid points to get "j".
As a result, taking a descrete derivative aproximation of "b" is **not** 
expected to *exactly* reproduce "j".

# Plots 

First lets look at "b1", which does not include the dipole field and so we
can isolate effects from that. It's divergence should be (close to) 0.

###LUHMANN1979
> ![](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled/divergence_B1.png)
> ![](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled2/divergence_B1.png)
> for luhmann1979, normalized divergence of B1 is almost always well within 10^-2 of the regions

###DIPTSUR2
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled/divergence_B1.png)
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled2/divergence_B1.png)
> for Diptsur2, normalized divergence of B1 is much higher, from 1 to 100%, across all distances.

Now lets put back the dipole field and look at the divergence of "b".

###LUHMANN1979
!> [](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled/divergence_B.png)
!> [](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled2/divergence_B.png)
> for luhmann1979, normalized divergence of B starts at arround 10^-2 then rapidly drops to near zero with distance,
> until very large distances where larger value appear in the normalized divergence, but not the relative divergence.
> Since we know the divergence of "b1" is so close to zero, this is almost entirely due to the dipole field.
> The large distances larger values of the normalized divergence are likely due to all the dipole partial derivatives being so close to zero.

###DIPTSUR2
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled/divergence_B.png)
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled2/divergence_B.png)
> for diptsur2, normalized divergence of B apears much the same as for B1

Now lets look at the error in amperes law.

###LUHMANN1979
> ![](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled/curlB1_and_J_percent_error.png)
> ![](LUHMANN1979/using_vtk/derivatives/20000101T000000/native_random_sampled2/curlB1_and_J_percent_error.png)
> The error starts out small at small distances, it gets noticably larger at around 5 Re (expected) but not that large,
> and then it seems to level out for very large distances. The change at 5 R_e is expected because the grid size gets larger,
> leading to larger discretization errors. The grid continues to get larger farther away
> I don't know why, but for large distances there appears to be a sharp line at around 10% error. 

###DIPTSUR2
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled/curlB1_and_J_percent_error.png)
> ![](DIPTSUR2/using_vtk/derivatives/20190902T063000/native_random_sampled2/curlB1_and_J_percent_error.png)
> Again, based on what I know, this computation of discrete curl B1 should reproduce J exactly. But the errors are actually
> much larger than even the discretization errors in the LUHMANN1979 one.
