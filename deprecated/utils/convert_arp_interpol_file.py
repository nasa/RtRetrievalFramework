#!/usr/bin/env python

import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

ils_cycle = 5.9 
num_footprint = 8

colors = range(1,1017)
num_colors = len(colors)
bands = ['O2', 'WCO2', 'SCO2']

input_file_format = '%s/ils_%s_Footprint%d.txt'

if (len(sys.argv) < 3):
    print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_directory> <output_file_base> [<in_nm>]\n"
    sys.exit(1)
    
input_directory  = sys.argv[1]
output_file_base = sys.argv[2]

# Label for file columns
file_label = ['PARAMETER', 'COEFFICIENT' ]
for band_idx in range(len(bands)):
    file_label.append('ILS_%d' % (band_idx+1))

for foot_num in range(1,num_footprint+1):
    print 'Processing footprint %d' % foot_num
    data_files = []
    for band_file in [ input_file_format % (input_directory, band_name, foot_num) for band_name in bands ]:
        print 'Reading %s' % os.path.basename(band_file)
        data_files.append(OCO_Matrix(band_file))
        
    file_rows = [ file_obj.dims[0] for file_obj in data_files ]

    nwvl = file_rows[0] / num_colors
    num_out_rows = 2 + num_colors + nwvl*(num_colors+1)

    ils_mat_data = numpy.zeros((num_out_rows, len(bands)+2), dtype=float)
    
    for row_idx in range(ils_mat_data.shape[0]):
        ils_mat_data[row_idx, 0] = row_idx + 1
        ils_mat_data[row_idx, 1] = 1
        
    for band_idx in range(len(bands)):
        row_idx = 0

        # width used by L2 code for convultion size
        ils_mat_data[row_idx, band_idx+2] = -data_files[band_idx].data[0,1] / ils_cycle
        row_idx += 1

        ils_mat_data[row_idx, band_idx+2] = num_colors
        row_idx += 1

        for color_idx in range(num_colors):
            ils_mat_data[row_idx, band_idx+2] = colors[color_idx]
            row_idx += 1

        for wvl_idx in range(nwvl):
            ils_mat_data[row_idx, band_idx+2] = data_files[band_idx].data[wvl_idx, 1]
            row_idx += 1


        for color_idx in range(num_colors):
            for wvl_idx in range(nwvl):
                data_index = nwvl*color_idx+wvl_idx
                #print 'index = ', data_index
                ils_mat_data[row_idx, band_idx+2] = math.log( data_files[band_idx].data[data_index, 2] )
                row_idx += 1

        
    last_dot = output_file_base.rfind('.')
    out_filename = '%s_%d%s' % (output_file_base[0:last_dot], foot_num, output_file_base[last_dot:])

    print 'Writing: %s' % out_filename
    file_obj = OCO_Matrix()    
    file_obj.data = ils_mat_data
    file_obj.file_id = 'Table ILS Data for Footprint %d' % (foot_num)
    file_obj.labels = file_label
    file_obj.write(out_filename, auto_size_cols=False, verbose=True)

