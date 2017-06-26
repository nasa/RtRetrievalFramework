import sys
import glob
import os

from full_physics_swig._swig_wrap import *

for i in glob.glob(os.path.dirname(__file__) + "/*.py"):
    mname = os.path.basename(i).split('.')[0]
    exec('from .%s import *' % mname)
