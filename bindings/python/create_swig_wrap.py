# This is a short script for generating the top level wrapper code that
# initializes everything for SWIG. You don't normally run this directly,
# it is used by python.am for automatically creating this file.

import os
import sys

print sys.argv[0]
tmpl_dir = os.path.dirname(sys.argv[0]) + "/"
swig_wrap_template = open(tmpl_dir + "swig_wrap.tmpl").read()
prototypes = "\n".join(["  void init_%s(void);" % i for i in sys.argv[2:]])
initcmd = "\n".join(["  init_extension_module(package, \"_%s\", init_%s);" % (i, i) for i in sys.argv[2:]])
with open(sys.argv[1], 'w') as swig_wrap_fo:
    swig_wrap_fo.write(swig_wrap_template.format(prototypes=prototypes,
                                                 initcmd=initcmd))


