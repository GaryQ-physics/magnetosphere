#sunspot_run.py

from config import conf
import util
import cutplane_plot_demo5 as cpd5
import cutplane_plot_demo2_animate as cpd2a
import earth_surf_dB.py as es


util.generate_TESTANALYTIC_files()
cpd5.main('TESTANALYTIC')

cpd2a.main('TESTANALYTIC')

es.tofile('TESTANALYTIC', (2000,1,1,1,1,0), para=True)
