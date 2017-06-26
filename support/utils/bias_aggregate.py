#! /usr/bin/env python
# Script to aggregate the output

from full_physics import *
import scipy.io
import numpy as np
version = "May 31, 2013"
usage = '''Usage:
   bias_aggregate [options] <output_file> <input_file> <input_file>...
   bias_aggregate -h | --help
   bias_aggregate -v | --version

This generates an aggregate matplotlib output file.

Options:
  -h --help         
     Print this message

  -v --version      
     Print program version
'''

args = docopt_simple(usage, version=version)
res = {}
for f in args.input_file:
    d = scipy.io.loadmat(f)
    for k in list(d.keys()):
        # Skip keys starting with "__"
        if(not k[0:2] == "__"):
            # Reshape data to have a extra dimension of size 1 in front
            s = [1]
            s.extend(d[k].shape)
            t = d[k].reshape(s)
            if(k in res):
                res[k] = np.append(res[k], t, axis = 0)
            else:
                res[k] = t
scipy.io.savemat(args.output_file, res)

