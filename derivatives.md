# Derivatives

The [total derivative](https://en.wikipedia.org/wiki/Total_derivative) of a vector field at a given point is a linear operator on the vectors at that point. This can be seen by arranging all the partial derivatives into [a matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant). An operator has a natural notion of "magnitude", given by the [operator norm](https://en.wikipedia.org/wiki/Operator_norm). In what follows, we will use 'NormD of F' as shorthand for the operator norm of the total derivative of F. Given a vector value function which in it's given units has a NormD well above machine precision, then other first derivative combinations (like divergence and curl) can be compared in a unit independent way to NormD. For example, if we have a field F where we want the divergence of F to be 0, then so long as NormD is well above machine precision, it makes sense to require div(F)/NormD(F)<<1.

For the SWMF output magnetic field and current field on the DIPTSUR2 run at 06:30, we will plot their divergence using the files original units for the field and R_E for the derivative. In addition, we will normalize it by using NormD and just the field itself.

The partial derivative were computed by symmetric difference. For points within 5 R_E of earth's center, a step size of 1/16 was used, i.e. points for the difference quotient were taken to be +1/16 and -1/16 away from the given point along each coordinate axis, for a total interval of 1/8. 
For points farther than 5 R_E from earth's center, a step size of 1/8 R_E was used.

# At time 20190902T063000

## J (MHD Current Field)
![](images/DIPTSUR2/20190902T063000/divergence_J.png)

>for all three scatter plots, each data point is associated to a point in space on the native BATS-R-US grid (chosen at random), and the x-axis is the distance of that point from the center of the earth.
>
>the y-axis differs for the 3 plots:
>
> - div(J) in units of (muA/m^2)/R_E
> - div(J)/norm(J) in units of 1/R_E
> - div(J)/NormD(J) which is unitless

## B1 (MHD Magnetic Field)
![](images/DIPTSUR2/20190902T063000/divergence_B1.png)

>for all three scatter plots, each data point is associated to a point in space on the native BATS-R-US grid (the same points as in the J case), and the x-axis is the distance of that point from the center of the earth.
>
>the y-axis differs for the 3 plots:
>
> - div(B1) in units of nT/R_E
> - div(B1)/norm(B1) in units of 1/R_E
> - div(B1)/NormD(B1) which is unitless

Now to check to the consistency of amperes law, since we expect J to be nonzero, we can compare curl(B) to mu0*J in a unit independent way without using NormD.
Both quantities are vector quantities, so we could either check all components seperately or use a vector difference. We will examin the vector difference first:

## Amperes law percent error
![](images/DIPTSUR2/20190902T063000/ampere_percent_error.png)

>The vertical line at rCurrents
>
>Each data point is associated to a point in space on the native BATS-R-US grid (the same points as before)
>
>Percent Error is given by 100.* norm(curl(B) - mu0*J) / norm(mu0*J)  and is plotted on the y-axis
>
>plotted on the x-axis is the distance of the corresponding point from the center of the earth in R_E.

the previous plots are using the field values at 6:30, when the magnetopause is closest to earth.
the following. For comparison, we now show data from 4:10.

# At time 20190902T041000

![](images/DIPTSUR2/20190902T041000/divergence_J.png)

![](images/DIPTSUR2/20190902T041000/divergence_B1.png)

![](images/DIPTSUR2/20190902T041000/ampere_percent_error.png)
