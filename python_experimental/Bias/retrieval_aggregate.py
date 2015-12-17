from full_physics import *
import shelve
import scipy.io

version = "August 30, 2013"
usage = '''Usage:
  retrieval_aggregate.py [options] <sounding_id>
  retrieval_aggregate.py -h | --help
  retrieval_aggregate.py -v | --version

Aggregate together the results.

Options:
  -h --help         
     Print this message

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
x_true = []
xco2_true = []
x_ret = []
xco2_ret = []
for fname in glob.glob(str(args.sounding_id) + "/retrieval_second*.shlv"):
    d = shelve.open(fname, "r")
    x_a = d["x_a"]
    cov_a = d["cov_a"]
    x_true.append(d["x_true"])
    xco2_true.append(d["xco2_true"])
    x_ret.append(d["x_retrieved"])
    xco2_ret.append(d["xco2_retrieved"])

res = {}
res["X_a"] = x_a
res["cov_a"] = cov_a
res["X_true"] = np.vstack(x_true)
res["XCO2_true"] = np.vstack(xco2_true)
res["X_retrieved"] = np.vstack(x_ret)
res["XCO2_retrieved"] = np.vstack(xco2_ret)
scipy.io.savemat(str(args.sounding_id) + "/l2_simulate.mat", res)
