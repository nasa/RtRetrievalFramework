from full_physics import *
import shelve
import cPickle

version = "August 29, 2013"
usage = '''Usage:
  create_draw_cov.py [options] <sounding_id> <output> <output_pkl>
  create_draw_cov.py -h | --help
  create_draw_cov.py -v | --version

This program looks through the retrievals that we could actually perform,
and from those determine the x_a and covariance of the real distribution.

Options:
  -h --help         
     Print this message

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)

x_true = []
for f in glob.glob(str(args.sounding_id) + "/draw_*.shlv"):
    din = shelve.open(f, "r")
    x_true.append(din["x_true"])
x_true = np.column_stack(x_true)
res = shelve.open(args.output)
res["x_a"] = np.mean(x_true, axis = 1)
res["cov_a"] = np.cov(x_true)
with open(args.output_pkl, "w") as f:
    cPickle.dump(np.cov(x_true), f)
