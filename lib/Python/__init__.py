import sys
import glob
import os

from try_swig_load import *

# Make sure we can safely import matplotlib without getting an error
# (see this module for details on this)
import safe_matplotlib_import

# Normal __all__ doesn't work here, because we also need to include the
# full_physics_swig stuff (with __all__, *only* the listed modules are
# included). So we just explicitly import the modules we want here.

for i in glob.glob(os.path.dirname(__file__) + "/*.py"):
    exec 'from ' + os.path.basename(i).split('.')[0] + ' import *'
