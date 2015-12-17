#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser
import numpy
import bisect
import contextlib
import copy

import GOSAT_File

from OCO_Matrix import OCO_Matrix

from matplotlib import pyplot as plt

FILE_ID_TMPL = 'GOSAT Spectrum for observation: %s'
FILE_LABELS = ['Radiance', 'Error', 'Wavelength', 'Wavenumber' ]
FILE_UNITS  = ['W/cm^2/sr/cm-1', 'W/cm^2/sr/cm-1', 'micron', 'cm-1']

BAND_RANGES = [ [12870, 13250], [5750, 6450], [4750, 5150] ]
BAND_ERROR_PIXELS = [ range(0,1000) + range(5565,6565),
                      range(0,1000) + range(7080,8080),
                      range(0,1000) + range(5565,6565) ]

HEADER_FILL_VALUES = {
    'spacecraft_alt'    : 0.0,
    'spacecraft_lat'    : 0.0,
    'spacecraft_lon'    : 0.0,
    'sounding_altitude' : 0.0,
    'relative_velocity' : 0.0,
    'polarization_angle': 0.0,
    }

def translate_header_name(curr_name):
    if curr_name in GOSAT_File.GEOMETRY_DATASETS.keys():
        return 'sounding_' + curr_name
    else:
        return curr_name

def write_spectra_files(output_dir, ids, data):
    for curr_id, curr_data in zip(ids, data):
        out_filename = '%s/%s.dat' % (output_dir, curr_id)
        used_data = []
        
        out_obj = OCO_Matrix()
        out_obj.labels = FILE_LABELS
        out_obj.units  = FILE_UNITS
        out_obj.file_id = FILE_ID_TMPL % (curr_id)

        out_obj.data = numpy.zeros((curr_data['Num_Points'], len(FILE_LABELS)), dtype=float)
        used_data.append('Num_Points')
        for column_idx, column_name in enumerate(FILE_LABELS):
            if curr_data.has_key(column_name):
                out_obj.data[:, column_idx] = curr_data[column_name]
                used_data.append(column_name)

        header_dict = copy.copy(curr_data)
        header_dict.update(HEADER_FILL_VALUES)
        for data_name, data_values in header_dict.items():
            if data_name in used_data:
                continue

            out_name = translate_header_name(data_name)

            try:
                if len(data_values.shape) == 0:
                    out_obj.header[out_name] = str(data_values)
                else:
                    out_obj.header[out_name] = ' '.join([ str(item) for item in iter(data_values)])
            except:
                try:
                    out_obj.header[out_name] = str(data_values)
                except:
                    print >>sys.stderr, 'Barfed on parsing %s with values: %s for header of file %s' % (data_name, data_values, out_filename)

        print 'Writing: %s' % out_filename
        out_obj.write(out_filename)

def filter_data_range(wavenumbers, spectrum):
    
    for curr_idx, curr_range in enumerate(BAND_RANGES):
        band_beg = float('inf')
        band_end = float('-inf')
        for pol_idx in range(len(GOSAT_File.POL_ORDER)):
            band_beg = min(band_beg, bisect.bisect_left(wavenumbers[curr_idx][:, pol_idx], curr_range[0]))
            band_end = max(band_beg, bisect.bisect_right(wavenumbers[curr_idx][:, pol_idx], curr_range[1]))

        wavenumbers[curr_idx] = wavenumbers[curr_idx][band_beg:band_end, :]
        spectrum[curr_idx]    = spectrum[curr_idx][band_beg:band_end, :]

def extract_spectra(hdf_file, output_dir, sens_corr=True, combine_pol=False, filter=False):

    print 'Opening HDF file: %s' % hdf_file
    with contextlib.closing(GOSAT_File.L1B(hdf_file, sensitivity_correction=sens_corr)) as l1b_obj:

        for obs_index in range(l1b_obj.get_num_observations()):
            obs_id = l1b_obj.get_observation_id(obs_index)

            print 'Reading data for observation: %s' % obs_id
            spectrum    = l1b_obj.get_swir_spectrum(obs_index)
            wavenumbers = l1b_obj.get_swir_wavenumbers(obs_index)
            geometry    = l1b_obj.get_swir_geometry(obs_index)
            time_stamp  = l1b_obj.get_time_stamp(obs_index)

            ids    = []
            data   = []

            if filter:
                filter_data_range(wavenumbers, spectrum)

            num_wns = [ band_wn.shape[0] for band_wn in wavenumbers ]

            wavelength_1d = numpy.zeros(numpy.sum(num_wns), dtype=float)
            wavenumber_1d = numpy.zeros(numpy.sum(num_wns), dtype=float)
            radiance_1d   = numpy.zeros(numpy.sum(num_wns), dtype=float)
            error_1d      = numpy.zeros(numpy.sum(num_wns), dtype=float)
            
            for pol_idx, pol_name in enumerate(GOSAT_File.POL_ORDER):
                start_pixels = numpy.zeros(len(spectrum), dtype=int)
                data_index = 0

                if not combine_pol:
                    radiance_1d[:] = 0.0
                    
                for band_idx, curr_wns, curr_spec in zip(range(len(spectrum)), wavenumbers, spectrum):
                    start_pixels[band_idx] = data_index + 1
                    pixel_index = 1
                    for inp_index in range(len(curr_spec)):
                        wavelength_1d[data_index] =  1e4/curr_wns[inp_index, pol_idx]
                        wavenumber_1d[data_index] =  curr_wns[inp_index, pol_idx]
                        radiance_1d[data_index]  += curr_spec[inp_index, pol_idx]
                        error_1d[data_index]      = numpy.std( radiance_1d[BAND_ERROR_PIXELS[band_idx]] )

                        pixel_index += 1
                        data_index  += 1

                if combine_pol and pol_idx == 0:
                    ids.append( '%s' % (obs_id) )
                else:
                    ids.append( '%s%s' % (obs_id, pol_name) )

                # Get average of spectra before setting
                if combine_pol and pol_idx == 1:
                    radiance_1d[:] /= len(GOSAT_File.POL_ORDER)
                
                if not combine_pol or pol_idx == 1:
                    
                    data_dict = { 'Num_Points': numpy.sum(num_wns),
                                  'Wavelength': wavelength_1d,
                                  'Wavenumber': wavenumber_1d,
                                  'Radiance'  : radiance_1d,
                                  'Error'     : error_1d,
                                  'Start_Pixels': start_pixels,
                                  'frame_time_stamp': time_stamp,
                                  }
                    for geom_name, geom_data in geometry.items():
                        if len(geom_data.shape) == 2:
                            data_dict[geom_name] = geometry[geom_name][:, pol_idx]
                        else:
                            data_dict[geom_name] = geometry[geom_name]

                    data.append( data_dict )

            write_spectra_files(output_dir, ids, data)

def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] hdf_filename [output_dir]")

    parser.add_option( "--raw", dest="sens_corr",
                       default=True,
                       action="store_false",
                       help="do not apply sensitivity correction to data"
                       )

    parser.add_option( "-c", dest="combine_pol",
                       default=False,
                       action="store_true",
                       help="combine polarization channels for each measurement/band"
                       )

    parser.add_option( "-f", dest="filter",
                       default=False,
                       action="store_true",
                       help="do not include out of band data in file"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if (len(args) < 1):
        parser.error("Need input_flename")

    hdf_file  = args[0]

    if not os.path.exists(hdf_file):
        parser.error("%s does not exist" % hdf_dump_file)

    if len(args) >= 2:
        output_dir = args[1]
    else:
        output_dir = '.'

    extract_spectra(hdf_file, output_dir, sens_corr=options.sens_corr, filter=options.filter, combine_pol=options.combine_pol)

if __name__ == "__main__":
    standalone_main()

