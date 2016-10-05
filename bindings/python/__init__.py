import sys
import glob
import os

from full_physics_swig._swig_wrap import *

for i in glob.glob(os.path.dirname(__file__) + "/*.py"):
    exec 'from ' + os.path.basename(i).split('.')[0] + ' import *'


