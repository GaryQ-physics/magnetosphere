#### configuration file to be imported by timeseries.py  ####
#### this is to be modified   ####

run = 'UNALT_DIPOLE4'

do_stats_summary = True
do_biotsavart_integral = True
do_coulomb_integral = True
do_probing = True
do_cutplane = True

integral_points = ['colaba','origin', (0.71875,0.09375,-3.71875), (-146.,-14.,-14.), (-156.,-124.,-124.)]
probe_points    = ['colaba','origin', (0.71875,0.09375,-3.71875), (-146.,-14.,-14.), (-156.,-124.,-124.)]

#  uncomment to overide rcut value from the default ( rCurrents of the corresponding run )
#rcut = 2.
