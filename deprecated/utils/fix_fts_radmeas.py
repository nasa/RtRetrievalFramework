#!/usr/bin/env python

import os
import math
import sys

from types import ListType
from OCO_Matrix import OCO_Matrix

def fix_fts_radmeas(input_file, output_file, scale, which_band=None):

    # Load existing file
    matrix_obj = OCO_Matrix(input_file)

    radiance_col = matrix_obj.labels_lower.index("radiance")

    # Swap wavelength and wavenumber column headers and units
    # since are reversed in FTS file output
    column1 = matrix_obj.labels_lower.index("wavelength")
    column2 = matrix_obj.labels_lower.index("wavenumber")

    labels_tmp = matrix_obj.labels[column2]
    units_tmp   = matrix_obj.units[column2]

    matrix_obj.labels[column2] = matrix_obj.labels[column1]
    matrix_obj.units[column2]  = matrix_obj.units[column1]
    
    matrix_obj.labels[column1] = labels_tmp
    matrix_obj.units[column1]  = units_tmp

    # Fix start pixels for individual band files
    if which_band != None:
        print "Fixing band pixel indexes"
        matrix_obj.pixels = [-2, -2, -2]

        matrix_obj.pixels[int(which_band)] = 0
        if(int(which_band) < 2):
            matrix_obj.pixels[int(which_band)+1] = matrix_obj.dims[0]-1

    # Scale radiance by specified amount
    if scale != None:
        print "Scaling radiances by %s" % scale
        for row_idx in range(matrix_obj.dims[0]):
            matrix_obj.data[row_idx][radiance_col] = matrix_obj.data[row_idx][radiance_col] * float(scale)

    print "Writing file: %s" % output_file
    matrix_obj.write(output_file,  auto_size_cols=False, verbose=True)
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    if type(scriptOptions) is ListType:
        fix_fts_radmeas(fileObj.filename, fileObj.filename, scriptOptions[0], scriptOptions[1])
    else:
        fix_fts_radmeas(fileObj.filename, fileObj.filename, scriptOptions)


def standalone_main():
    if (len(sys.argv) < 5):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<input_matrix_file> <output_output_file> <scale>\n"
        sys.exit(1)

    input_matrix_file  = sys.argv[1]
    output_matrix_file = sys.argv[2]
    scale              = sys.argv[3]
    which_band         = sys.argv[4]

    fix_fts_radmeas(input_matrix_file, output_matrix_file, scale, which_band)

if __name__ == "__main__":
    standalone_main()
