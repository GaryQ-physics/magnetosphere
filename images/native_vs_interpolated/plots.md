# Intro

These plots are for the DIPTSUR2 run at 6:30.
The 'old interpolated way' is using the VTK library interpolation, not Kameleon.

b and b1 are in nanotesla, j in microamperes per meter squared. Divergence uses units of 1/R_E

all points in the native grid are used except those at the edges/boundaries of a block,
and those within 1.8 R_E from earths center.

"jR" here is the attempt to reconstruct "j" from "b1". This is done by 
taking curl_b1, and dividing by mu_0 = 8014.956840
(vacuum permeability in the apropriate units)

In all the plots that follow, the dots are color coded based on
whether the point is on the day side or the night side.
Night side is Blue, Day side is Orange.

## divergence of b

![](div_b.png)

## norm (magnitude) of b1

![](norm_b1.png)

## divergence of b1

![](div_b1.png)

## norm (magnitude) of j

![](norm_j.png)

## divergence of j

![](div_j.png)

## norm (magnitude) of jR

![](norm_jR.png)

## error in jR

![](jR_error.png)

> this is the norm (magnitude) of the vector difference between jR and j,
> and thus is still in microamperes per meter squared.

# Closer Look

It appears from the above "error in jR" and "norm (magnitude) of jR" plots
that the error is tends to be several orders of magnitude smaller than the magnitude.

In order to quantify this, we will plot the ratio of the two:

# jR fractional error

![](jR_fractional_error.png)

For completeness, also plot the ratio of divergence of b1 and magnitude of b1:

# relative divergence of B1

![](rel_div_b1.png)



