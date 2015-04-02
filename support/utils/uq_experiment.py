#!/usr/bin/env python

import os
import sys

import h5py
import numpy

from full_physics import *

class UqExperiment(object):

    def __init__(self, config_filename, output_file):
        # Load Lua state, and config object
        self.ls, self.lua_config = l2_lua.load_lua_config(config_filename)

        # Where results are saved
        self.output_file = output_file

        # Now create everything
        self.lua_config.do_config(self.lua_config)

    def run(self):
        l2_lua.run_retrieval(self.ls, self.output_file)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="")

    parser.add_argument("lua_config", help="Lua configuration file")
    parser.add_argument("output_file", help="Filename for L2 output")

    # Parse command line arguments
    args = parser.parse_args()

    uq = UqExperiment(args.lua_config, args.output_file)
    uq.run()
