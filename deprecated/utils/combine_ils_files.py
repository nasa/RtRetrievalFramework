#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

# Format for file number
band_names = ['O2', 'WCO2', 'SCO2']
inp_file_format = 'ils_%s_Footprint%d.txt'
out_file_format = 'ils_table_%d.dat'
num_footprints = 8
num_pixels = 1016
num_wvl = 200
num_rows = num_pixels*num_wvl

header_additionals = {
    'function_type':      'TABLE',
    'num_ils_parameters': str(num_pixels),
    'num_ils_wndepend':   str(num_wvl),
    'interpolation':      '100 100 100',
                       }

file_labels = [ 'ILS_PIXELS' ]
for band_idx in range(len(band_names)):
    file_labels.append('ILS_DELTA_LAMBDA_%d' % (band_idx+1))
    file_labels.append('ILS_RESPONSE_%d' % (band_idx+1))

parser = OptionParser(usage="usage: %prog [options] <input_dir> <output_dir>")

# Parse command line arguments
(options, args) = parser.parse_args()

if (len(args) < 1):
    parser.error("Need input dir at least")
input_dir = args[0]

if len(args) >= 2:
    output_dir = args[1]
else:
    output_dir = '.'

for fp_idx in range(num_footprints):
    output_data_filename = ('%s/' + out_file_format) % (output_dir, fp_idx+1)
    print 'Processing footprint: %d' % (fp_idx+1)

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'Instrument Line Shape parameters for frame %d' % (fp_idx+1)

    file_data = numpy.zeros((num_rows, len(band_names)*2 + 1), dtype=float) 

    for band_idx in range(len(band_names)):
        input_data_filename = ('%s/' + inp_file_format) % (input_dir, band_names[band_idx], fp_idx+1)
        print 'Reading from: %s' % input_data_filename
        in_mat_obj = OCO_Matrix(input_data_filename)

        if num_rows != in_mat_obj.dims[0]:
            raise IOError('Num rows expected: %d does not match number read: %d from file: %s'  (num_rows, in_mat_obj.dims[0], input_data_filename))

        for row_idx in range(num_rows):
            if band_idx == 0:
                file_data[row_idx, 0] = in_mat_obj.data[row_idx, 0]

            for col_num in range(1,3):
                src_idx = col_num
                dst_idx = col_num+(band_idx*2)
                if row_idx == 0:
                    print 'Set value %d -> %d' % (src_idx, dst_idx)
                file_data[row_idx, dst_idx] = in_mat_obj.data[row_idx, src_idx]

    out_mat_obj.labels = file_labels
    out_mat_obj.data = file_data

    for (key_name, key_val) in header_additionals.iteritems():
        out_mat_obj.header[key_name] = key_val

    print 'Writing to: %s' % output_data_filename
    out_mat_obj.write(output_data_filename)
