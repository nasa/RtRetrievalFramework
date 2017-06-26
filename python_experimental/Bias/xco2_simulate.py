from full_physics import *
import cPickle

version = "August 19, 2013"

usage = '''Usage:
  xco2_simulate [options] <dist_cov> <config_file> <sounding_id> <output>
  xco2_simulate -h | --help
  xco2_simulate -v | --version

This program is used to empirical calculate the XCO2 bias. This does a 
single run of a simulated Level 2 to return a XCO2 retrieved value. We
start with the x_a supplied by an existing Level 2 run. We then perturb
this randomly according to a supplied distribution covariance matrix. We
run the forward model and add random measurement error according to s_e
read from the Level 2 run. Finally we do a Level 2 retrieval.

Note that the distribution covariance matrix and the apriori covariance
matrix are not in general the same matrix.

If there is a successful retrieval, the output is a pickled result that
we can then collect together over a number of runs.

Options:
  -h --help         
     Print this message

  --lua-config=s
     Specify the Lua config file to use. Default is normal 
     gosat/config/config.lua [default: ]

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
lua_config = args.lua_config
if(lua_config == ""):
    lua_config = None

res = {}
# Start with a Level 2 run
l2run = L2Run.create_apriori_state(args.config_file, 
                                   sounding_id = str(args.sounding_id), 
                                   sounding_id_index = 0, 
                                   lua_config = lua_config)
# Gives us x_a and the apriori covariance matrix
x_a = np.copy(l2run.state_vector.state)
l2run_cov_a = np.copy(l2run.state_vector.state_covariance)
with open(args.dist_cov) as f:
    dist_cov = cPickle.load(f)
res["XCO2_a"] = l2run.xco2
res["X_a"] = np.copy(l2run.state_vector.state)
print "XCO2 at x_a:", l2run.xco2
print "Apriori state vector:"
print l2run.state_vector

# Generate a random "true" state

# We have diag(v.T * cov_a * v) = w, so
# Gives random distribution with covariance cov_a, mean x_a
w, v = np.linalg.eigh(dist_cov)
x = x_a + np.dot(v, np.random.randn(w.shape[0]) * np.sqrt(w))
l2run.state_vector.update_state(x, l2run_cov_a)
res["X_true"] = np.copy(l2run.state_vector.state)
res["XCO2_true"] = l2run.xco2
print "XCO2 at X_true:", l2run.xco2
print "'True' state vector:"
print l2run.state_vector

# Simulate Level 1b radiance data, adding in random noise
rad = []
for i in range(3):
    sr = l2run.forward_model.radiance(i, True).spectral_range
    l1rad = l2run.forward_model.level_1b.radiance(i)
    sr = sr.convert(l1rad.units)
    noise = np.random.randn(sr.data.shape[0]) * l1rad.uncertainty
    t = sr.data + noise
    rad.append(SpectralRange(sr.data + noise, sr.units, l1rad.uncertainty))

# Update what the solver is looking at to match the simulated L1b 
# radiance data
l2run.forward_model.level_1b = Level1bCache(l2run.forward_model.level_1b)
for i in range(3):
    l2run.forward_model.level_1b.set_radiance(i, rad[i])

# Retrieve
l2run.solver.solve(x_a, x_a, l2run_cov_a)
l2run.forward_model.state_vector.update_state(l2run.solver.x_solution,
                                              l2run.solver.aposteriori_covariance)
res["X_retrieved"] = np.copy(l2run.state_vector.state)
res["XCO2_retrieved"] = l2run.xco2
print "Retrieved XCO2:", l2run.xco2
print "X_retrieved"
print l2run.state_vector
with open(args.output, "w") as f:
    cPickle.dump(res, f)

