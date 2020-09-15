#sunspot_run.py

from config import conf
import util
import cutplane_plot_demo2b as cpd2b
#import cutplane_plot_demo2_animate as cpd2a
import earth_surf_dB as es


#util.generate_TESTANALYTIC_files()
cpd2b.main('DIPTSUR2')

#cpd2a.main('DI')

es.tofile('DIPTSUR2', (2019,9,2,6,30,0,0), para=True)
