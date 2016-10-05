#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix
from H5Dump_Util import H5Dump_Util

parser = OptionParser(usage="usage: %prog [options] hdf_filename [output_dir]")

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


h5dump_obj = H5Dump_Util(hdf_file, debug=False)

(num_bands, num_soundings, num_pixels, num_coefs) = h5dump_obj.Get_Dims('/InstrumentHeader/snr_coef')

print 'Reading: %s' % hdf_file
snr_data = h5dump_obj.Get_Dataset_Values('/InstrumentHeader/snr_coef')

for snd_idx in range(num_soundings):
    output_filename = '%s/snr_coefs_%d.dat' % (output_dir, snd_idx+1)

    file_data = numpy.zeros((num_pixels, num_bands*num_coefs+1), dtype=float)

    for pix_idx in range(num_pixels):
        column_labels = []
        col_idx = 0
        
        column_labels.append('COLOR')
        file_data[pix_idx, col_idx] = pix_idx+1
        col_idx += 1

        for band_idx in range(num_bands):
            column_labels.append('C_PHOTON_%d' % (band_idx+1))
            file_data[pix_idx, col_idx] = snr_data[band_idx, snd_idx, pix_idx, 0]
            col_idx += 1

        for band_idx in range(num_bands):
            column_labels.append('C_BACKGROUND_%d' % (band_idx+1))
            file_data[pix_idx, col_idx] = snr_data[band_idx, snd_idx, pix_idx, 1]
            col_idx += 1

    matrix_obj = OCO_Matrix()
    matrix_obj.data = file_data
    matrix_obj.labels = column_labels
    matrix_obj.file_id = 'SNR Coefficients Extracted from %s' % hdf_file

    print 'Writing coefs to: %s' % output_filename
    matrix_obj.write(output_filename)
    
    
