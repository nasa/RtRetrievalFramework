#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import os
import sys
import logging
from argparse import ArgumentParser

import h5py
import numpy as np

from full_physics import *
import full_physics

# Some constants that normally would be obtained from Lua config or static input
NUM_PIXEL = 1016
ILS_HALF_WIDTH = [ DoubleWithUnit(4.09e-04, "um"), 
                   DoubleWithUnit(1.08e-03, "um"),
                   DoubleWithUnit(1.40e-03, "um") ]
ILS_INTERPOLATE = False

DESC_BAND_NAMES = ["A-Band", "WC-Band", "SC-Band"]
HDF_BAND_NAMES = ["o2", "weak_co2", "strong_co2"]

# This was the latest NADIR at the time script was made, use this if no L1B is specified
DEFAULT_L1B_FILE = "/oco2/product/Ops_B4300_r01/2014/09/14/L1bSc/oco2_L1bScND_01082a_140914_B4300_140915132020.h5"

DEFAULT_SOUNDING = 3

def convolve_data(l1b_obj, high_res_wl, high_res_rad, band_idx, sounding_pos=DEFAULT_SOUNDING, in_wavenumber=False):
    # Load ILS table, make sure we are passing double sized arrays
    delta_lambda = np.array(l1b_obj['/InstrumentHeader/ils_delta_lambda'][:][band_idx, sounding_pos, :, :], dtype=np.float64)
    response = np.array(l1b_obj['/InstrumentHeader/ils_relative_response'][:][band_idx, sounding_pos, :, :], dtype=np.float64)
    disp_coeff = np.array(l1b_obj['/InstrumentHeader/dispersion_coef_samp'][:][band_idx, sounding_pos, :], dtype=np.float64)

    # Divide out and flip the order if input file defined in wavenumbers 
    if in_wavenumber:
        high_res_wl = old_div(1e4,high_res_wl[::-1])
        high_res_rad = high_res_rad[::-1]

    # Create dispersion object
    disp_l2 = DispersionPolynomial(disp_coeff, np.zeros(disp_coeff.shape, dtype=bool), "micron", DESC_BAND_NAMES[band_idx], NUM_PIXEL, True)
    wavelength = disp_l2.pixel_grid.data

    # Make sure range of wavelength grids will result in real data
    logging.debug("High Resolution Range = [%s, %s]" % (high_res_wl[0], high_res_wl[-1]))
    logging.debug("Dispersion Range = [%s, %s]" % (wavelength[0], wavelength[-1]))

    if not (high_res_wl[0] < wavelength[0] and high_res_wl[-1] > wavelength[-1]):
        logging.warning('High resolution wavelength range: [%s, %s] must be larger than dispersion wavelength range: [%s, %s]' % (high_res_wl[0], high_res_wl[-1], wavelength[0], wavelength[-1]))

    # Perform convolution
    ils_table = IlsTableLog(wavelength, delta_lambda, response, DESC_BAND_NAMES[band_idx], HDF_BAND_NAMES[band_idx], ILS_INTERPOLATE)
    ils_conv = IlsConvolution(disp_l2, ils_table, ILS_HALF_WIDTH[band_idx])
    conv_rad = ils_conv.apply_ils(high_res_wl, high_res_rad, list(range(0,NUM_PIXEL)))

    # Convert to ph / s / micron
    conv_rad = conv_rad / disp_l2.pixel_grid.photon_to_radiance_factor().value

    return wavelength, conv_rad

def convolve_data_from_file(l1b_file, rad_file, band_idx, sounding_pos=DEFAULT_SOUNDING, in_wavenumber=False):
    rad_data = np.loadtxt(rad_file)
       
    # Load high resolution radiance, make sure we use double sized data
    high_res_wl = np.array(rad_data[:, 0], dtype=np.float64)
    high_res_rad = np.array(rad_data[:, 1], dtype=np.float64) 

    return convolve_data(l1b_file, high_res_wl, high_res_rad, sounding_pos, in_wavenumber)

def display_result(conv_wl, conv_rad):
    for wl, rad in zip(conv_wl, conv_rad):
        print("%.06e %.06e" % (wl, rad))

if __name__ == "__main__":
    parser = ArgumentParser(description='Convolves a simple text file with the L2 ils table')

    parser.add_argument('rad_file', metavar='FILENAME',
            help='Radiance file with 2 columns, one with spectral grid and one with radiance data')
    parser.add_argument('band_index', type=int,
            help='Which band index, [0,2]')
    parser.add_argument('--sounding_pos', '-s', type=int, default=DEFAULT_SOUNDING,
            help='Index of the sounding posistion to use for ILS convolution')
    parser.add_argument('--l1b_file', '-l', metavar='FILENAME', default=DEFAULT_L1B_FILE,
            help='Alternative L1B file to get ILS from other than: %s' % DEFAULT_L1B_FILE)
    parser.add_argument('--wavenumbers', default=False, action='store_true',
            help='Specify this flag if the grid in rad_file is defined using wavenumbers')
    parser.add_argument('--verbose', '-v', default=False, action='store_true',
            help='Output verbose debugging information')

    args = parser.parse_args()

    if args.verbose:
        log_util.init_logging(logging.DEBUG)
    else:
        log_util.init_logging(logging.INFO)

    # Debugging to let us make sure we are using the correct L2 library
    logging.debug("L2 Library in use: %s" % full_physics.__file__)

    l1b_obj = h5py.File(args.l1b_file, 'r')

    conv_wl, conv_rad = convolve_data(l1b_obj, args.rad_file, args.band_index, args.sounding_pos, in_wavenumber=args.wavenumbers)
    display_result(conv_wl, conv_rad)
