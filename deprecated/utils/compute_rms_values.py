#!/usr/bin/env python

import os
import sys
import re
from time import gmtime, strftime

import numpy

from OCO_Matrix import OCO_Matrix
import OCO_Statistics

if (len(sys.argv) < 3):
    print "usage:\n\t", sys.argv[0], "<dir_pairs_file> <output_data_filename> [rms_compare_file]\n"
    print "Outputs the RMS value for each pair of testcase directories in the input file"
    sys.exit(1)

tcDirPairsFilename = sys.argv[1]
outputDataFilename = sys.argv[2]

if len(sys.argv) >= 4:
    rmsCompareFile = sys.argv[3]
else:
    rmsCompareFile = 'rad_conv.dat'

pairsFObj = open(tcDirPairsFilename, 'r')

pairsList = pairsFObj.readlines()

pairsFObj.close()

rms_data = numpy.zeros((len(pairsList) * 3, 2), dtype=float)
pair_index = 0
for pair in pairsList:
    pair =  pair.replace('#', '').strip()

    pair_parts = pair.split()

    if len(pair_parts) > 1:
        (runDir, refDir) = pair_parts[0:2]
    else:
        runDir = pair_parts[0]
        refDir = os.path.dirname(runDir) + '/std_output'
    
    runCompareFilename = "%s/out/%s" % (runDir, rmsCompareFile)
    refCompareFilename = "%s/out/%s" % (refDir, rmsCompareFile)

    if not (os.path.exists(runCompareFilename) and os.path.exists(refCompareFilename)):
        print "%s not found in run and reference dirs for %s" % (rmsCompareFile, runDir)
        continue

    print runCompareFilename, refCompareFilename

    runFileObj = OCO_Matrix(runCompareFilename)
    refFileObj = OCO_Matrix(refCompareFilename)

    if (rmsCompareFile.find("high_res") == -1):
        windows = [0, 1024, 2048, runFileObj.dims[0]]
    else:
        windows = [0, runFileObj.dims[0]]


    index = runFileObj.labels_lower.index("radiance")

    # Plot each window separately
    for w in range(0, len(windows) - 1):

        rundata = runFileObj.data[windows[w]:windows[w+1], index]
        refdata = refFileObj.data[windows[w]:windows[w+1], index]
        residual, rms = OCO_Statistics.RMS(rundata, refdata)

        print '%s %d %f' % (os.path.basename(runDir), w, rms)
        
        rms_data[pair_index][0] = w
        rms_data[pair_index][1] = rms
     
        pair_index = pair_index + 1


out_mat_obj = OCO_Matrix()
out_mat_obj.file_id = 'Relative RMS Difference for: %s' % rmsCompareFile
out_mat_obj.labels = ['Window_Num', 'RMS_Diff']
out_mat_obj.data = rms_data

out_mat_obj.write(outputDataFilename)
