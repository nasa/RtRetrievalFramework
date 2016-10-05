from full_physics import *
import shelve

# Script to make a draw from a random distribution, and also to generate the
# Level 1 data for it.

version = "August 29, 2013"
usage = '''Usage:
  create_draw.py [options] <sounding_id> <output>
  create_draw.py -h | --help
  create_draw.py -v | --version

This program does a draw for the given sounding, and calculates the
Level 1 from the forward model.

Options:
  -h --help         
     Print this message

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
din = shelve.open(str(args.sounding_id) + "/cov_initial.shlv", "r")
x_a = din["x_a"]
cov_a = din["cov_a"]
cov_a_original = din["cov_a_original"]

basedir = "/groups/algorithm/l2_fp/oco2_simulator_test/automated/tags/B3.5.00_aerosol_testing_0/Orbit022/"

l2run = L2Run.create_from_existing_run(basedir + "oco_Orbit022.config", 
                                       str(args.sounding_id), 0,
                                       "config_quicklook_aerosol_test_0.lua")
res = shelve.open(args.output)

# Generate a random "true" state

# We have diag(v.T * cov_a * v) = w, so
# Gives random distribution with covariance cov_a, mean x_a
w, v = np.linalg.eigh(cov_a)
x = x_a + np.dot(v, np.random.randn(w.shape[0]) * np.sqrt(w))
l2run.state_vector.update_state(x, cov_a_original)
res["x_true"] = np.copy(l2run.state_vector.state)
res["xco2_true"] = l2run.xco2

# Simulate Level 1b radiance data, adding in random noise
rad = []
scale = 1.0
for i in range(3):
    sr = l2run.forward_model.radiance(i, True).spectral_range
    l1rad = l2run.forward_model.measured_radiance(i).spectral_range
    sr = sr.convert(l1rad.units)
    noise = np.random.randn(sr.data.shape[0]) * l1rad.uncertainty * scale
    t = sr.data + noise
    rad.append(SpectralRange(sr.data + noise, sr.units, 
                             l1rad.uncertainty * scale))

res["radiance_simulated"] = rad
res.close()

