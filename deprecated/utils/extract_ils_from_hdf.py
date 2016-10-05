#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser
import numpy
import contextlib

import h5py

from OCO_Matrix import OCO_Matrix

def extract_ils_from_hdf(hdf_file, output_dir, reverse=False):
    print 'Opening HDF file: %s' % hdf_file
    with contextlib.closing(h5py.File(hdf_file, 'r')) as hdf_obj:

        reported_soundings = []
        for snd_idx, rep_val in enumerate(hdf_obj['Metadata']['ReportedSoundings']):
            if rep_val > 0:
                reported_soundings.append( snd_idx + 1 )

        ils_delta_lambda      = hdf_obj['InstrumentHeader']['ils_delta_lambda']
        ils_relative_response = hdf_obj['InstrumentHeader']['ils_relative_response']

        num_bands          = ils_delta_lambda.shape[0]
        num_ils_parameters = ils_delta_lambda.shape[2]
        num_ils_wndepend   = ils_delta_lambda.shape[3]

        labels = [ 'ILS_PIXELS' ]
        for band_num in range(1, num_bands+1):
            labels.append('ILS_DELTA_LAMBDA_%d' % band_num)
            labels.append('ILS_RESPONSE_%d' % band_num)

        for snd_idx, sounding_id in enumerate(reported_soundings):
            output_filename = os.path.join(output_dir, 'ils_%d.dat' % sounding_id)
            print 'Extracting data for sounding %d into %s' % (sounding_id, output_filename)

            ils_file = OCO_Matrix()
            ils_file.data = numpy.zeros((num_ils_parameters * num_ils_wndepend, num_bands*2+1), dtype=float)

            row_beg = 0
            for color_index in range(num_ils_parameters):
                print 'Extracting color %d' % (color_index+1)
                row_end = row_beg+num_ils_wndepend
                ils_file.data[row_beg:row_end, 0] = color_index + 1

                for band_idx, col_idx in zip(range(num_bands), range(1,1+2*num_bands,2)):
                    ils_file.data[row_beg:row_end, col_idx]   = ils_delta_lambda[band_idx, snd_idx, color_index, :]

                    if reverse:
                        ils_file.data[row_beg:row_end, col_idx+1] = ils_relative_response[band_idx, snd_idx, color_index, :][::-1]
                    else:
                        ils_file.data[row_beg:row_end, col_idx+1] = ils_relative_response[band_idx, snd_idx, color_index, :]

                row_beg = row_end


            ils_file.file_id = 'Instrument Line Shape parameters for sounding posistion %d' % sounding_id

            ils_file.labels = labels

            ils_file.header['function_type'] = 'TABLE'
            ils_file.header['interpolation'] = '100 100 100'
            ils_file.header['num_ils_parameters'] = num_ils_parameters
            ils_file.header['num_ils_wndepend']   = num_ils_wndepend

            print 'Writing to %s' % output_filename
            ils_file.write(output_filename, verbose=True)


def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] hdf_filename [output_dir]")

    parser.add_option( "-r", "--reverse_response", dest="reverse_response",
                       default=False,
                       action="store_true",
                       help="Reverse response values, for use by badly written HDF files")

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

    extract_ils_from_hdf(hdf_file, output_dir, options.reverse_response)

if __name__ == "__main__":
    standalone_main()
