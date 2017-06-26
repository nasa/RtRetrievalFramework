#!/usr/bin/env python
# This run a Level 2 retrieval, but using python as the top level. 
# In particular, we can use a python config file.

import os
import imp
import re
from full_physics import *

def run_python_config(config_filename, output_file):
    """Runs a python configurations which is expected to load a Lua configuration itself
       and create a global in it's name space named "ls" 
    """
    config = imp.load_source('config', config_filename)
    run_retrieval(config.ls, output_file)


def run_lua_config(config_filename, output_file):
    "Runs a Lua configuration as is"
    
    ls, lua_config = l2_lua.load_lua_config(config_filename)
    lua_config.do_config(lua_config)
    l2_lua.run_retrieval(ls, output_file)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
    usage="""usage: %prog <config file> <output file>

    This run a Level 2 retrieval, but using python as the top level. 
    In particular, we use a python config file, which can use python
    classes to prototype changes to L2.
    """)

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) == 2:
        # Fully expand filenames so we can chdir without worrying
        # about not using paths as supplied on command line
        config_filename, output_file = [ os.path.realpath(fn) for fn in args ]
        
        if re.search('\.py$', config_filename):
            run_python_config(config_filename, output_file)
        elif re.search('\.lua$', config_filename):
            run_lua_config(config_filename, output_file)
        else:
            raise ValueError("Configuration filename has unrecognized extension, not .py or .lua: %s" % config_filename)
    else:
        parser.error("Need to specify all the arguments")
