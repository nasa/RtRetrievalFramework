#!/usr/bin/env python

from builtins import str
from builtins import zip
from builtins import range

import os
import sys
import logging
from argparse import ArgumentParser

from glob import glob
from datetime import datetime

import h5py
import numpy as np

from full_physics import log_util
from full_physics.l2_input import L2InputFile
from full_physics import HdfFile, SpectralWindowRange, Level1bFts, DispersionPolynomial

sys.path.append(os.path.dirname(sys.argv[0]))
from offline_ils import convolve_data, DEFAULT_SOUNDING

DEFAULT_FTS_CONFIG_FILE = os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]), '../../input/fts/config/config.lua'))

FTS_STATIC_INPUT_FILE = os.path.realpath(os.path.join(os.path.dirname(sys.argv[0]), '../../input/fts/input/l2_fts_static_input.h5'))

# Arbitrarily picked a nadir from the the last day of data for B7 
REFERENCE_L1B_FILE = "/data/oco2/ops/product/Ops_B7000_r02/2015/11/09/L1bSc/oco2_L1bScND_07218a_151109_B7000r_151220114932.h5"

DATASETS_TO_COPY_FROM_L1B = ['/InstrumentHeader/ils_delta_lambda', '/InstrumentHeader/ils_relative_response',
        '/InstrumentHeader/snr_coef', '/InstrumentHeader/dispersion_coef_samp']

BAND_NAMES = ['o2', 'weak_co2', 'strong_co2']

logger = logging.getLogger()

def spectral_window(static_input_file=FTS_STATIC_INPUT_FILE):
    static_inp = HdfFile(static_input_file)
    win_ranges = static_inp.read_double_with_unit_3d("Spectral_Window/microwindow")
    return SpectralWindowRange(win_ranges)

def create_sounding_id(obs_id, l1b_fts, sounding_pos=DEFAULT_SOUNDING):
    dt = datetime.utcfromtimestamp(l1b_fts.time(0).unix_time)
    sounding_id = '%s%03d%d' % (dt.strftime('%Y%m%d%H%M%S'), obs_id, (sounding_pos+1))
    return sounding_id

def create_datasets(obs_ids, output_file, l1b_file):
    "Create the datasets for the output file sized according to the number of frames being processed, will only have one sounding id per frame"

    # Copy instrument related data from the reference L1B file
    for ds_name in DATASETS_TO_COPY_FROM_L1B:
        src_ds = l1b_file[ds_name]
        output_file.create_dataset(ds_name, dtype=src_ds.dtype, data=src_ds[:])

    num_frames = len(obs_ids)
    output_file.create_dataset('/FootprintGeometry/footprint_altitude', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_azimuth', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_zenith', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_latitude', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_longitude', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_solar_azimuth', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_solar_zenith', (num_frames,1,3), dtype=float)
    output_file.create_dataset('/FootprintGeometry/footprint_time_tai93', (num_frames,1,3), dtype=float)

    rv = output_file.create_dataset('/FrameGeometry/relative_velocity', (num_frames,), dtype=float)
    rv[:] = 0

    output_file.create_dataset('/FrameHeader/frame_time', (num_frames,), dtype=float)
    output_file.create_dataset('/FrameHeader/frame_time_stamp', (num_frames,), dtype='S25')

    output_file.create_dataset('/SoundingGeometry/sounding_id', (num_frames,1), dtype=np.int64)

    lw = output_file.create_dataset('/SoundingGeometry/sounding_land_water_indicator', (num_frames,1), dtype=int)
    lw[:,:] = 0

    for band_name in BAND_NAMES:
        output_file.create_dataset('/SoundingMeasurements/radiance_%s' % band_name, (num_frames,1,1016), dtype=float)
        output_file.create_dataset('/SoundingMeasurements/wavelength_%s' % band_name, (num_frames,1,1016), dtype=float)

def store_l1b_data(sounding_id, snd_idx, l1b_fts, output_file, l1b_file, sounding_pos=DEFAULT_SOUNDING):

    for band_idx in range(len(BAND_NAMES)):
        fp_time = l1b_fts.time(band_idx)

        output_file['/FootprintGeometry/footprint_altitude'][snd_idx, 0, band_idx] = l1b_fts.altitude(band_idx).value
        output_file['/FootprintGeometry/footprint_azimuth'][snd_idx, 0, band_idx] = l1b_fts.sounding_azimuth(band_idx).value
        output_file['/FootprintGeometry/footprint_zenith'][snd_idx, 0, band_idx] = l1b_fts.sounding_zenith(band_idx).value
        output_file['/FootprintGeometry/footprint_latitude'][snd_idx, 0, band_idx] = l1b_fts.latitude(band_idx).value
        output_file['/FootprintGeometry/footprint_longitude'][snd_idx, 0, band_idx] = l1b_fts.longitude(band_idx).value
        output_file['/FootprintGeometry/footprint_solar_azimuth'][snd_idx, 0, band_idx] = l1b_fts.solar_azimuth(band_idx).value
        output_file['/FootprintGeometry/footprint_solar_zenith'][snd_idx, 0, band_idx] = l1b_fts.solar_zenith(band_idx).value
        output_file['/FootprintGeometry/footprint_time_tai93'][snd_idx, 0, band_idx] = fp_time.pgs_time

        output_file['/SoundingGeometry/sounding_id'][snd_idx, 0] = int(sounding_id)

        # Measured radiance from FTS insturment to convolve
        high_res_rad = l1b_fts.radiance(band_idx).data

        # Come up with grid associated with radiance data from dispersion coefficients
        disp_coeff = l1b_fts.spectral_coefficient(band_idx)
        disp_flag = np.zeros((disp_coeff.rows), dtype=bool)
        number_pixel = high_res_rad.shape[0]
        disp_poly = DispersionPolynomial(disp_coeff.value, disp_flag, disp_coeff.units, BAND_NAMES[band_idx], number_pixel, False)
        high_res_wl = disp_poly.pixel_grid.data

        conv_wl, conv_rad = convolve_data(l1b_file, high_res_wl, high_res_rad, band_idx, in_wavenumber=True, sounding_pos=sounding_pos)

        output_file['/SoundingMeasurements/radiance_%s' % BAND_NAMES[band_idx]][snd_idx, 0, :] = conv_rad
        output_file['/SoundingMeasurements/wavelength_%s' % BAND_NAMES[band_idx]][snd_idx, 0, :] = conv_wl

    output_file['/FrameHeader/frame_time'][snd_idx] = fp_time.pgs_time
    output_file['/FrameHeader/frame_time_stamp'][snd_idx] = bytes(fp_time)

def convert_fts_data(fts_run_dir, output_filename, l1b_filename=REFERENCE_L1B_FILE, config_filename=DEFAULT_FTS_CONFIG_FILE,
                     window_filename=DEFAULT_FTS_CONFIG_FILE, sounding_pos=DEFAULT_SOUNDING):
    # Open files needed for processing
    output_file = h5py.File(output_filename, 'w')
    l1b_file = h5py.File(l1b_filename, 'r')

    # Get per sounding information from files in the populated FTS run directory
    obs_ids = [ int(sid) for sid in np.loadtxt(os.path.join(fts_run_dir, 'sounding_id.list')) ]
    
    spec_a_files = [l.strip() for l in open(os.path.join(fts_run_dir, 'spectrum_a.list'), 'r').readlines()]
    spec_b_files = [l.strip() for l in open(os.path.join(fts_run_dir, 'spectrum_b.list'), 'r').readlines()]
    # Initialize output file datasets
    create_datasets(obs_ids, output_file, l1b_file)

    # Load the run directory configuration file to determine location of run log
    config_fn = glob(os.path.join(fts_run_dir, '*.config'))[0]
    config_file = L2InputFile(config_fn)
    config_sec = config_file.rootNode.get_section('input->InputProductFiles')[0]
    runlog_file = config_sec.get_keyword_value('RunlogFile')

    logger.debug('Run Log = %s' % runlog_file)

    spec_win = spectral_window(window_filename)

    # Load the FTS L1B object for each sounding, writing that information to the output file
    for idx, (obs_id, spec_a, spec_b) in enumerate(zip(obs_ids, spec_a_files, spec_b_files)):
        spectra_names = [spec_b, spec_a, spec_a]

        l1b_fts = Level1bFts(runlog_file, spectra_names, spec_win.spectral_bound)
        sounding_id = create_sounding_id(obs_id, l1b_fts,sounding_pos)

        logger.debug('Obs ID = %s, Sounding ID = %s' % (obs_id, sounding_id))
        logger.debug('Spectrum A = %s' % spec_a)
        logger.debug('Spectrum B = %s' % spec_b)

        store_l1b_data(sounding_id, idx, l1b_fts, output_file, l1b_file, sounding_pos=sounding_pos)

if __name__ == "__main__":
    parser = ArgumentParser(description='Converts FTS data into ')

    parser.add_argument('fts_dir', metavar='FILENAME')

    parser.add_argument('--output_file', '-o', metavar='FILENAME', required=True,
        help='File to write convolved FTS data in OCO format')

    parser.add_argument('--l1b_file', '-l', metavar='FILENAME', default=REFERENCE_L1B_FILE,
        help='L1B file to get ILS information from')

    parser.add_argument('--window_file', '-w', metavar='FILENAME', default=FTS_STATIC_INPUT_FILE,
        help='static input file to get window information from')

    parser.add_argument('--config_file', '-c', metavar='FILENAME', default=DEFAULT_FTS_CONFIG_FILE,
        help='static input file to get window information from')

    parser.add_argument('--sounding_pos', '-s', metavar='INT', default=DEFAULT_SOUNDING, type=int, choices=range(8),
        help='specify which footprint ils to use (zero-based index)')

    parser.add_argument('--verbose', '-v', default=False, action='store_true',
        help='Output verbose debugging information')

    args = parser.parse_args()

    if args.verbose:
        log_util.init_logging(logging.DEBUG)
    else:
        log_util.init_logging(logging.INFO)

    convert_fts_data(args.fts_dir, args.output_file, args.l1b_file, args.config_file, args.window_file, sounding_pos=args.sounding_pos)
