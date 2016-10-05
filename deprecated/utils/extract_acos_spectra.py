#!/usr/bin/env python

import os
import sys
import time
import contextlib
from optparse import OptionParser

import numpy
import h5py

import ACOS_File
from OCO_Matrix import OCO_Matrix

FILE_ID_TMPL = 'ACOS Spectrum for observation: %s'
FILE_LABELS = ['Radiance', 'Noise', 'Wavelength', 'Wavenumber' ]
FILE_UNITS  = ['W/cm^2/sr/cm-1', 'W/cm^2/sr/cm-1', 'micron', 'cm-1']

FILE_COL_INDEXES = {}
for col_index, label_name in enumerate(FILE_LABELS):
    FILE_COL_INDEXES[label_name] = col_index

def write_ascii_from_hdf(acos_l1b_obj, sounding_info_dict, ascii_obj, sounding_id, average_name=None):
    ascii_obj.file_id = '%s spectra for sounding id: %s' % (acos_l1b_obj.instrument_name.upper(), sounding_id)
    ascii_obj.labels = FILE_LABELS
    ascii_obj.units = FILE_UNITS
   
    ascii_obj.header = sounding_info_dict

    ascii_obj.header['sounding_id'] = sounding_id
    ascii_obj.header['hdf_l1b_filename'] = os.path.realpath(acos_l1b_obj.filename)

    radiance_data = acos_l1b_obj.get_radiance_data(sounding_id, average=average_name)
    wavenumbers   = acos_l1b_obj.get_wavenumbers(sounding_id, average=average_name)
    wavelengths   = acos_l1b_obj.get_wavelengths(sounding_id, average=average_name)
    error_data    = acos_l1b_obj.get_error_data(sounding_id, average=average_name)

    data_len = acos_l1b_obj.get_channel_counts(sounding_id, average=average_name)
    num_rows = numpy.sum(data_len)

    ascii_obj.data = numpy.zeros((num_rows, len(FILE_LABELS)), dtype=float)

    start_index = 0
    ascii_obj.pixels = []
    for band_idx, band_len in enumerate(data_len):
        ascii_obj.pixels.append(start_index)

        ascii_obj.data[start_index:start_index+band_len, FILE_COL_INDEXES['Radiance']] = radiance_data[band_idx][:]
        ascii_obj.data[start_index:start_index+band_len, FILE_COL_INDEXES['Noise']] = error_data[band_idx][:]
        ascii_obj.data[start_index:start_index+band_len, FILE_COL_INDEXES['Wavenumber']] = wavenumbers[band_idx][:]
        ascii_obj.data[start_index:start_index+band_len, FILE_COL_INDEXES['Wavelength']] = wavelengths[band_idx][:]
        start_index += band_len
    

def extract_acos_spectra(acos_l1b_file, output_dir, sounding_id_list=None, average_pol=False):

    if average_pol:
        average_name = 'Polarization'
    else:
        average_name = None
    
    with contextlib.closing(ACOS_File.L1B(acos_l1b_file)) as acos_l1b_obj:
        if sounding_id_list == None or len(sounding_id_list) == 0:
            sounding_id_list = acos_l1b_obj.get_sounding_ids(add_polarization=not average_pol)

        for sounding_id in numpy.ravel(sounding_id_list):
            ascii_obj = OCO_Matrix()

            sounding_info_dict = acos_l1b_obj.get_sounding_info_dict(sounding_id, ignore_missing=True, average=average_name)

            write_ascii_from_hdf(acos_l1b_obj, sounding_info_dict, ascii_obj, sounding_id, average_name)

            output_filename = os.path.join(output_dir, '%s.dat' % sounding_id)
            print 'Writing %s' % output_filename
            ascii_obj.write(output_filename, default_precision=16)

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] hdf_filename [output_dir]")
   
    parser.add_option( "-s", "--sounding_id", dest="sounding_id_list",
                       metavar="INT",
                       action="append",
                       help="desired sounding id for extraction")

    parser.add_option( "-a", "--avergage_pol", dest="average_pol",
                       action="store_true",
                       default=False,
                       help="average polarizations")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('Need to at least specify output file')

    acos_l1b_file = args[0]

    if len(args) < 2:
        output_dir = './'
    else:
        output_dir = args[1]

    extract_acos_spectra(acos_l1b_file, output_dir, options.sounding_id_list, average_pol=options.average_pol)

if __name__ == "__main__":
    standalone_main()
