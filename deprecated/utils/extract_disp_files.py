#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

# Format for file number
grep_cmd = 'grep -A 28 DispersionCoefSamp'
disp_sect_beg = 'DATASET "DispersionCoefSamp"'
disp_sect_skip = 3
disp_sect_end = '}'

num_coefs = 3
num_bands = 3

file_labels = ['PARAMETER', 'DISP_1', 'DISP_2', 'DISP_3' ]

parser = OptionParser(usage="usage: %prog [options] hdf_dump_file [output_dir]")

# Parse command line arguments
(options, args) = parser.parse_args()

if (len(args) < 1):
    parser.error("Need input_flename")

hdf_dump_file  = args[0]

if not os.path.exists(hdf_dump_file):
    parser.error("%s does not exist" % hdf_dump_file)

if len(args) >= 2:
    output_dir = args[1]
else:
    output_dir = '.'

print 'Opening: %s' % hdf_dump_file
dump_obj = os.popen('%s %s' % (grep_cmd, hdf_dump_file), 'r')

found_disp_sect = False
finished_file = False
dispersion_coefs = [ [] for bidx in range(num_bands) ]
while not finished_file:
    line_str = dump_obj.readline().strip()

    if found_disp_sect:

        if line_str.find(disp_sect_end) >= 0:
            print 'End of dispersion section reached'
            finished_file = True
        else:
            (indicator, values_str) = line_str.split(':')
            indexes = indicator.replace('(', '').replace(')', '').split(',')
            coefs = values_str.split(',')[0:num_coefs]
            band_idx = int(indexes[0])
            dispersion_coefs[band_idx].append(coefs)
    else:
        if line_str.find(disp_sect_beg) >= 0:
            print 'Found dispersion section'
            for skip_idx in range(disp_sect_skip):
                dump_obj.readline()
            found_disp_sect = True
dump_obj.close()

if finished_file and not found_disp_sect:
    raise IOError('Never found dispersion section in %s' % hdf_dump_file)

num_disp_files = len(dispersion_coefs[0])
for disp_index in range(num_disp_files):
    output_data_filename = '%s/dispersion_%d.dat' % (output_dir, disp_index+1)

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'Dispersion parameters for frame %d' % (disp_index+1)

    file_data = numpy.zeros((num_coefs, num_bands + 1), dtype=float)

    for coef_idx in range(num_coefs):
        file_data[coef_idx][0] = coef_idx + 1       

    for band_idx in range(1, num_bands+1):
        for coef_idx in range(num_coefs):
            file_data[coef_idx][band_idx] = float( dispersion_coefs[band_idx-1][disp_index][coef_idx] )

    out_mat_obj.labels = file_labels
    out_mat_obj.data = file_data

    print 'Writing: %s' % output_data_filename
    out_mat_obj.write(output_data_filename)
