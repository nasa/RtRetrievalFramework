#!/usr/bin/env python

from types import ListType
import re
import time
import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

def create_mean_psurf(runlog_file, psurf_file):

    print 'runlog_file = ', runlog_file
    print 'psurf_file = ', psurf_file
    
    runlog_fobj = open(runlog_file, "r")

    header_cols = runlog_fobj.readline().split()

    pout_col = header_cols.index('pout')

    pouts = []
    for runlog_line in runlog_fobj.readlines():
        runlog_parts = runlog_line.split()
        pouts.append(float(runlog_parts[pout_col]))
        
    runlog_fobj.close()
    
    avg_psurf = mean(pouts) * 1e2

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = "Mean surface pressure from runlog file: %s" % runlog_file
    out_mat_obj.labels = ['LEVEL', 'PSURF']
    out_mat_obj.data = ones((1, 2), dtype=float)
    out_mat_obj.data[0, 1] = avg_psurf
    out_mat_obj.write(psurf_file)


def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    create_mean_psurf(fileObj.filename, scriptOptions)

def standalone_main():
    if (len(sys.argv) < 3):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<runlog_file> <out_psurf_file>\n"
        sys.exit(1)

    runlog_filename   = sys.argv[1]
    psurf_filename    = sys.argv[2]

    create_mean_psurf(runlog_filename, psurf_filename)

if __name__ == "__main__":
    standalone_main()
