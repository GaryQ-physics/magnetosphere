#### configuration file to be imported by timeseries.py  ####
#### this is to be modified   ####

run = 'UNALT_DIPOLE'

do_stats_summary = False
do_cutplane = False

points = [
            {   'point' : "colaba",
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : False,
            },

#            {   'point' : "origin",
#                'do_biotsavart_integral' : True,
#                'do_coulomb_integral' : True,
#                'do_probing' : False,
#            },

#            {   'point' : (0.71875,0.09375,-3.71875),
#                'do_biotsavart_integral' : True,
#                'do_coulomb_integral' : True,
#                'do_probing' : True,
#                'probe_vars' : ['b1'],
#            },
#
#            {   'point' : (-146.,-14.,-14.),
#                'do_biotsavart_integral' : True,
#                'do_coulomb_integral' : True,
#                'do_probing' : True,
#                'probe_vars' : ['b1'],
#            },

            {   'point' : (-156.,-124.,-124.),
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : True,
                'probe_vars' : ['b1'],
            },
        ]


#  uncomment to overide rcut value from the default ( rCurrents of the corresponding run )
#rcut = 2.
