from full_physics import *
import numpy as np
from nose.tools import *
from matplotlib.pylab import *
import os

# This is a set of test data that we used when comparing against to
# OCO-2 simulator at CSU
bdir = os.path.dirname(__file__) + "/"
sdir = bdir + "simulator_result/"
sid_name = "2010090912004075"
met_name = "../input/oco2_sim_met.h5"
l1b_name = "../input/oco2_sim_l1b.h5"
scene_file = "../input/oco2_sim_scene.h5"
r = L2Run(bdir + "config/config_orbit_sim_match.lua", sid_name,met_name, 
          l1b_name, scene_file = scene_file)

# Read in the spectral domain that we will be comparing against. This is the
# first column in the text file (after skipping the header row)
sd = []
for b in range(3):
    sd.append(SpectralDomain(np.genfromtxt(sdir + "hires_10774_%d.txt" % (b+1), 
                                           skip_header = 1)[:,0], 
                             Unit("nanometer")))

def test_solar_spectrum():
    '''Compare OCO-2 simulator and JPL solar spectrum'''
    for b in range(3):
        jpl_spec = r.solar_model[b].solar_spectrum(sd[b])
        # Data is 5th column, after skipping 5 header rows
        csu_spec = np.genfromtxt(sdir + "solar_spectrum_10774_%d.txt" % (b + 1),
                                 skip_header = 5)[:,4]
        per_diff = (jpl_spec.value - csu_spec) / max(csu_spec) * 100.0
        assert max(abs(per_diff)) < 0.15

def test_plot_solar_spectrum():
    '''Plot out the difference. Not really a test'''
    for b in range(3):
        jpl_spec = r.solar_model[b].solar_spectrum(sd[b])
        # Data is 5th column, after skipping 5 header rows
        csu_spec = np.genfromtxt(sdir + "solar_spectrum_10774_%d.txt" % (b + 1),
                                 skip_header = 5)[:,4]
        per_diff = (jpl_spec.value - csu_spec) / max(csu_spec) * 100.0
        subplot(3,1,b+1)
        title("Difference full solar spectrum band %d" % (b+1))
        xlabel("Wavelength")
        ylabel("Difference %")
        plot(sd[b].wavelength(), per_diff)
    savefig("solar_diff.png")
    
