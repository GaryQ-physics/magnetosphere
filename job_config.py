#### configuration file to be imported by timeseries.py  ####
#### this is to be modified   ####

run = 'DIPTSUR2'

do_stats_summary = False
do_cutplane = False

points = [
            {   'point' : "colaba",
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : False,
                'do_gap_region' : True,
            },

#            {   'point' : "origin",
#                'do_biotsavart_integral' : True,
#                'do_coulomb_integral' : True,
#                'do_probing' : False,
#            },

            {   'point' : (0.71875,0.09375,-3.71875),
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : True,
                'probe_vars' : ['b1','b'],
                'do_gap_region' : True,

            },

            {   'point' : (-146.,-14.,-14.),
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : True,
                'probe_vars' : ['b1','b'],
                'do_gap_region' : True,
            },

            {   'point' : (-156.,-124.,-124.),
                'do_biotsavart_integral' : True,
                'do_coulomb_integral' : True,
                'do_probing' : True,
                'probe_vars' : ['b1','b'],
                'do_gap_region' : True,
            },
        ]


#  uncomment to overide rcut value from the default ( rCurrents of the corresponding run )
#rcut = 2.
