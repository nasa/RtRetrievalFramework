#!/usr/bin/env python

import os
import math
import sys
from optparse import OptionParser

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

# Format for file number
count_fmt = '-%02d'

parser = OptionParser(usage="usage: %prog [options] input_filename [output_basename]")

parser.add_option( "-m", "--max_files", dest="max_files",
                   metavar="NUM",
                   type="int",
                   help="maximum number of files to extract from input_filename"
                   )

# Parse command line arguments
(options, args) = parser.parse_args()

if (len(args) < 1):
    parser.error("Need input_flename")

input_fort_file  = args[0]

if not os.path.exists(input_fort_file):
    parser.error("%s does not exist" % input_fort_file)

if len(args) >= 2:
    output_base_fmt = args[1]
    
    ext_start = output_base_fmt.rfind('.')
    if ext_start > 0:
        ext_str = output_base_fmt[ext_start:]
        output_base_fmt = output_base_fmt.replace(ext_str, (count_fmt + ext_str))
    else:
        output_base_fmt += count_fmt
else:
    output_base_fmt = 'extracted-%s%s.dat' % (os.path.basename(input_fort_file), count_fmt)

print 'Opening: %s' % input_fort_file
fort_obj = open(input_fort_file, "r")

file_number = 1
file_line_num = 1
finished_file = False
while not finished_file and (options.max_files == None or file_number <= options.max_files):
    size_line = fort_obj.readline()
    size_parts = size_line.split()

    if len(size_parts) != 3:
        finished_file = True
        break

    file_line_num += 1

    num_rows = int(size_parts[0])
    num_cols = int(size_parts[1])
    num_head = int(size_parts[2])

    print 'Found section of size %d x %d' % (num_rows, num_cols)

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = "Extracted from %s" % input_fort_file

    header_count = 0
    for head_index in range(num_head):
        fort_head_line = fort_obj.readline()
        print 'Header %02d: %s' % fort_head_line
        out_mat_obj.header['Comment_%02d' % header_count] = '"%s"' % fort_head_line
        header_count += 1
        file_line_num += 1

    file_data = numpy.zeros((num_rows, num_cols), dtype=float)

    curr_row = 0
    curr_col = 0       
    for row_index in range(num_rows):
        data_line = fort_obj.readline().split()
        file_line_num += 1

        for data_val in data_line:
            if curr_col >= num_cols:
                curr_row += 1
                curr_col = 0

            file_data[curr_row, curr_col] = float(data_val)
            curr_col += 1

    out_mat_obj.data = file_data

    output_data_filename = output_base_fmt % file_number

    print 'Writing: %s' % output_data_filename
    out_mat_obj.write(output_data_filename)

    file_number += 1

fort_obj.close()

