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

        # Update the statevector with the values from the UQ file
        self.load_uq_sv()

    def load_uq_sv(self):
        uq_obj = HdfFile(self.lua_config.spectrum_file)
        frame_idx = self.lua_config.l1b_sid_list(self.lua_config).frame_number
        initial_sv = uq_obj.read_double_2d("/StateVector/sampled_state_vectors")[:, frame_idx]

        ig_orig = self.lua_config.initial_guess
        ig_v = fp.InitialGuessValue()
        ig_v.initial_guess = initial_sv
        ig_v.apriori = ig_orig.apriori
        ig_v.apriori_covariance = ig_orig.apriori_covariance
        ig_new = fp.CompositeInitialGuess()
        ig_new.add_builder(ig_v)
        self.ls.globals["initial_guess"] = ig_new
        self.lua_config.state_vector.update_state(initial_sv)


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
