#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import os
import sys
import itertools
from contextlib import closing

import h5py
import numpy

from full_physics import *

logger = FpLogger(FpLogger.INFO)
Logger.set_implementation(logger)

def log_info(info_str):
    logger.write(FpLogger.INFO, info_str)

def log_flush():
    logger.flush(FpLogger.INFO)

class StateVectorError(BaseException):
    pass

class UqExperiment(object):

    def __init__(self, config_filename, output_file):
        # Load Lua state, and config object
        self.ls, self.lua_config = l2_lua.load_lua_config(config_filename)

        # Seed the random number generator with the sounding id number
        # to ensure consistent results when rerun
        numpy.random.seed(int(self.lua_config.sid_string) % (2**32))

        # Where results are saved
        self.output_file = output_file

        # Now create everything
        self.lua_config.do_config(self.lua_config)

        # Number of spectrometers
        self.num_spec = self.lua_config.number_pixel.rows(self.lua_config.number_pixel)

    def update_initial_guess(self, new_sv):
        'Modify the inital guess to values different than those loaded by Lua'

        ig_orig = self.lua_config.initial_guess
        ig_v = fp.InitialGuessValue()
        ig_v.initial_guess = new_sv
        ig_v.apriori = ig_orig.apriori
        ig_v.apriori_covariance = ig_orig.apriori_covariance
        ig_new = fp.CompositeInitialGuess()
        ig_new.add_builder(ig_v)
        self.ls.globals["initial_guess"] = ig_new
        try:
            self.lua_config.state_vector.update_state(new_sv)
        except RuntimeError as exc:
            raise StateVectorError(exc.message)

    def calculate_uq_radiance(self):
        'Compute radiances for use by UQ L1B. Noise is added based on the L1B noise model'

        # Load the true state vector from UQ file for the Measured RT calculation
        orig_sv = self.lua_config.state_vector.state.copy()
        with closing(HdfFile(self.lua_config.spectrum_file)) as uq_obj:
            frame_idx = self.lua_config.l1b_sid_list(self.lua_config).frame_number
            true_sv = uq_obj.read_double_2d("/StateVector/sampled_state_vectors")[:, frame_idx]

            try:
                self.update_initial_guess(true_sv)
            except StateVectorError as exc:
                self.compare_svs(exc.message)
                raise exc

            if uq_obj.has_object("/SpectralParameters/sampled_radiance_offset"):
                all_rad_offset = uq_obj.read_double_2d("/SpectralParameters/sampled_radiance_offset")[:, frame_idx]
            else:
                all_rad_offset = None
                log_info("Could not find /SpectralParameters/sampled_radiance_offset in file: %s\n" % uq_obj.file_name)
                log_flush()
  
        # Change SpectralWindows so that all radiances points are computed
        spec_win_range = self.lua_config.spec_win.range_array
        orig_spec_win_value = spec_win_range.value.copy()
        new_spec_win_value = spec_win_range.value.copy()

        orig_bad_samp_mask = self.lua_config.spec_win.bad_sample_mask.copy()
        new_bad_samp_mask = self.lua_config.spec_win.bad_sample_mask.copy()

        for spec_idx in range(int(self.num_spec)):
            # Use all samples (pixels), + 1 because the range is not inclusive
            new_spec_win_value[spec_idx, 0, 0] = 0
            new_spec_win_value[spec_idx, 0, 1] = self.lua_config.number_pixel(spec_idx) + 1

        new_bad_samp_mask[:, :] = False

        spec_win_range.value = new_spec_win_value
        self.lua_config.spec_win.range_array = spec_win_range
        self.lua_config.spec_win.bad_sample_mask = new_bad_samp_mask

        # Calculate radiances for use by UQ object
        skip_jacobian = True
        for spec_idx in range(int(self.num_spec)):
            uq_spectrum = self.lua_config.fm.config.forward_model.radiance(spec_idx, skip_jacobian)
            uq_spec_range = uq_spectrum.spectral_range.convert(Unit("Ph sec^{-1} m^{-2} sr^{-1} um^{-1}"))

            if all_rad_offset is not None:
                pixel_range = self.lua_config.fm.config.forward_model.pixel_range(spec_idx)
                band_rad_offset = all_rad_offset[pixel_range.first():pixel_range.last()+1]
            else:
                band_rad_offset = numpy.zeros(uq_spec_range.data.shape, dtype=uq_spec_range.data.dtype)

            # Add noise to radiances using nose model uncertainty for gaussian, use
            # sounding id number as seed for reproducability
            uncertainty = self.lua_config.l1b.noise_model.uncertainty(spec_idx, uq_spec_range.data)
            noisified_radiance = numpy.random.normal(uq_spec_range.data + band_rad_offset, uncertainty)
            noisified_range = SpectralRange(noisified_radiance, uq_spec_range.units)

            # Save noisified spectra into the L1B class
            self.lua_config.l1b.set_radiance(spec_idx, noisified_range)

        # Reset spectral window range to original value
        spec_win_range.value = orig_spec_win_value
        self.lua_config.spec_win.range_array = spec_win_range
        self.lua_config.spec_win.bad_sample_mask = orig_bad_samp_mask

        # Reset initial guess to values loaded by Lua
        self.update_initial_guess(orig_sv)

    def compare_svs(self, err_message=None):
        if err_message != None:
            print(err_message, file=sys.stderr)

        file_sv_names = []
        with closing(h5py.File(self.lua_config.spectrum_file, "r")) as uq_obj:
            file_sv_names = uq_obj['/StateVector/state_vector_names'][:]
        config_sv_names = self.lua_config.state_vector.state_vector_name

        name_len = 0
        for name in itertools.chain(file_sv_names, config_sv_names):
            name_len = max(name_len, len(name))

        hdr_format = "    %"+str(name_len)+"s %"+str(name_len)+"s"
        line_format = "% 3d %"+str(name_len)+"s %"+str(name_len)+"s"

        print(hdr_format % ("Scence File SV Name", "Config SV Name"), file=sys.stderr)
        print("-" * (name_len*2+6), file=sys.stderr)
        for idx, (file_name, config_name) in enumerate(itertools.zip_longest(file_sv_names, config_sv_names)):
            print(line_format % (idx, file_name, config_name), file=sys.stderr)

    def run(self):
        # Calculate the radiances need for the retrieval
        self.calculate_uq_radiance()

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
