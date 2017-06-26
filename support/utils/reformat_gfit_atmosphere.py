#!/usr/bin/env python
from __future__ import print_function
from full_physics.oco_matrix import OcoMatrix

from types import ListType
import sys
import os
import re
import numpy
import glob

mav_col_extract = [ ['Pres', 1.01325e5],
                    ['Height', 1],
                    ['Temp', 1 ],
                    ['1h2o', 1 ],
                    ['1co2', 1 ], 
                    ['1o2',  1 ],
                    ['1o3',  1 ],
                    ['1ch4',  1 ],]

all_col_names = ['Pressure', 'Height', 'Temperature', 'H2O', 'CO2', 'O2', 'O3', 'CH4']
all_unit_names = ['Pa', 'km', 'K', 'VMR', 'VMR', 'VMR', 'VMR', 'VMR']


def reformat_gfit_atmosphere(mav_file, out_file, next_spec_srch=None):

    print('Reading mav data from %s' % mav_file)

    spec_line_start = 1

    mav_data = []
    mav_fobj = open(mav_file, "r")
    file_line_idx = 0

    if next_spec_srch == None:
        found_spec = True
    else:
        found_spec = False
    for mav_line in mav_fobj.readlines():
        line_parts = mav_line.split()
        mav_data.append( line_parts )

        if mav_line.find('Next Spectrum:') >= 0 and next_spec_srch != None:
            if re.search(next_spec_srch, mav_line):
                spec_line_start = file_line_idx
                found_spec = True
                
        file_line_idx += 1

    if not found_spec:
        raise ValueError('Could not find next spectrum search string: %s in mav file: %s' % (next_spec_srch, mav_file))

    print('Processing for', ' '.join(mav_data[spec_line_start]))

    mav_size_row   = spec_line_start + 1
    mav_header_row = mav_size_row + 4

    try:
        (num_skip, num_cols, num_rows) = [int(val) for val in mav_data[mav_size_row]]
    except:
        mav_header_row = 0
        num_skip = -2
        num_cols = len(mav_data[0])
        num_rows = len(mav_data)

    print() 

    print("Skip: %d, Cols %d, Rows: %d" % (num_skip, num_cols, num_rows))

    mav_beg_row = mav_size_row + num_skip + 2
    mav_end_row = mav_beg_row + num_rows - 3

    mav_all_cols = mav_data[mav_header_row]

    print("Column names:", mav_all_cols)

    out_col_idx = 0
    output_data_matrix = numpy.zeros((mav_end_row-mav_beg_row+1, len(all_col_names)), dtype=float)

    for (curr_mav_col, scale) in mav_col_extract:
        print('Processing column:', curr_mav_col)
        mav_col_idx = mav_all_cols.index(curr_mav_col)
        row_idx = mav_end_row-mav_beg_row

        for mav_row_data in mav_data[mav_beg_row:mav_end_row+1]:
            new_col_data = float(mav_row_data[mav_col_idx]) * float(scale)
            output_data_matrix[row_idx, out_col_idx] = output_data_matrix[row_idx, out_col_idx] + new_col_data
            row_idx -= 1

        out_col_idx += 1

    print('Writing output file %s' % out_file)
    out_mat_obj = OcoMatrix()
    out_mat_obj.file_id = 'GFIT Atmospheric State modified from: %s' % (mav_file)
    out_mat_obj.dims = [len(output_data_matrix), len(all_col_names)]
    out_mat_obj.labels = all_col_names
    out_mat_obj.units = all_unit_names
    out_mat_obj.data = output_data_matrix
    out_mat_obj.write(out_file)

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    opt_mav_filename = ''
    next_spec_srch = None
    if type(scriptOptions) is ListType:
        opt_mav_filename = scriptOptions[0]
        next_spec_srch   = scriptOptions[1]
    else:
        opt_mav_filename = scriptOptions

    print(opt_mav_filename)

    try:
        input_mav_file = (glob.glob( os.path.expanduser(opt_mav_filename) ))[0]
    except:
        return IOError('Could not find mav file: %s' % opt_mav_filename)
    
    reformat_gfit_atmosphere(input_mav_file, fileObj.filename, next_spec_srch)

def standalone_main():
    if len(sys.argv) < 3:
        print("usage:\n\t", os.path.basename(sys.argv[0]), "<gfit .mav file> <output_filename> [next_spec_search_re]\n")
        print("Creates a L2 atmosphere.dat file from a GFIT input files")
        sys.exit(1)

    mav_file = sys.argv[1]
    out_file = sys.argv[2]

    next_spec_srch = None
    if len(sys.argv) >= 4:
        next_spec_srch = sys.argv[3]
    
    reformat_gfit_atmosphere(mav_file, out_file, next_spec_srch)

if __name__ == "__main__":
    standalone_main()



