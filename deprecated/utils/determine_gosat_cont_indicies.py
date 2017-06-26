#!/usr/bin/env python

import os
import sys
import numpy
import bisect
import contextlib

import GOSAT_File


WN_RANGES = ( ( (12968.0, 12976.0), (13186.5, 13190.0) ),
              ( ( 6227.0,  6228.0), ( 6229.0,  6230.0) ),
              ( 4817.0, 4854.0, 4896.1 ) )

hdf_file = sys.argv[1]
with contextlib.closing(GOSAT_File.L1B(hdf_file, correct_disp=True)) as l1b_obj:
    wavenumbers = l1b_obj.get_swir_wavenumbers(1)

    for band_idx, band_ranges in enumerate(WN_RANGES):
        range_strings = []
        indv_indexes = []
        for curr_range in band_ranges:
            if type(curr_range) is tuple:
                (beg_wn, end_wn) = curr_range
                beg_idx = bisect.bisect(wavenumbers[band_idx][:, 0], beg_wn)
                end_idx = bisect.bisect(wavenumbers[band_idx][:, 0], end_wn)
                range_strings.append( 'range(%d, %d)' % (beg_idx, end_idx) )
            else:
                wn_idx = bisect.bisect(wavenumbers[band_idx][:, 0], curr_range)
                indv_indexes.append(wn_idx)

        if len(indv_indexes) > 0:
            print '%s' % indv_indexes,
        if len(indv_indexes) > 0 and len(range_strings)  > 0:
            print  ' + ',
        if len(range_strings):
            print ' + '.join(range_strings),
        print
