# Biot-Savart and MHD Comparison

A set of points were selected at random from the native SWMF grid. At each of these points, the following quantities were calculated and recorded:

- *B1_simulation* : the MHD magnetic field excluding the Earth dipole field ("b1" in SWMF output files)
- *B_biotsavart* : the magnetic field calculated due to the MhD current density ("j" in SWMF output files)

The points on the Day and Night side were plotted seperately in two different figures. The plots show percent error vs distance from Earth center. The percent error of each cartesian component is calculated with: 
- the corresponding component of *B1_simulation* used as the "expected value"
- the corresponding component of *B_biotsavart*  used as the "observed value"

All distances are in Earth radii, all magnetic fields are in nanoTesla.

## Nightside
![](images/DIPTSUR2/night_side.png) 
>3 scatterplots, one for each cartesian component. One dot for each sampled gridpoint on the nightside.
>
>datapoints color coded based on the absolute value of the corresponding component of *B1_simulation* (or *value*)
>
>small includes only points with *value* < 100 nT.
>mid includes only points with 100 nT <= *value* <= 1000 nT.
>large includes only points with 1000 nT <= *value*.

## Dayside
![](images/DIPTSUR2/day_side.png) 
>3 scatterplots, one for each cartesian component. One dot for each sampled gridpoint on the dayside.
>
>datapoints color coded based on the absolute value of the corresponding component of *B1_simulation* (or *value*)
>
>small includes only points with *value* < 100 nT.
>mid includes only points with 100 nT <= *value* <= 1000 nT.
>large includes only points with 1000 nT <= *value*.

Note: Points outside the MHD domain, because they are say too close to earth, have *B1_simulation* = 0. These are not included in the plots, as they would give values of +/- infinity.

## Summary of Results and Conclusion
Overall, the immediate take away is the large percent errors. They were getting to even an order of magnitude off (1000%) in some places for the mid values (>100nT). The large values (>1000nT) only occur less than 5 R\_E in. For the large values on the dayside, there is a noticeable trend downward from 100% to 10% or so as you go from 2 to 5 R\_E. For the nightside, there are too few large values to tell. However, since there doesn't appear to be a clear trend among the mid-range points which extend much farther than 5 R\_E, these results do not support the idea that the large percent errors are mainly due to a 'shielding effect' from the earth being a conductor.
