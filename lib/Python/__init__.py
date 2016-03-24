from __future__ import absolute_import
import sys
import glob
import os
import re

from .try_swig_load import *

# Make sure we can safely import matplotlib without getting an error
# (see this module for details on this)
from . import safe_matplotlib_import

# Normal __all__ doesn't work here, because we also need to include the
# full_physics_swig stuff (with __all__, *only* the listed modules are
# included). So we just explicitly import the modules we want here.

for i in glob.glob(os.path.dirname(__file__) + "/*.py"):
    mname = os.path.basename(i).split('.')[0]
    if(not re.match('.*_test', mname)):
        exec('from .%s import *' % mname)
