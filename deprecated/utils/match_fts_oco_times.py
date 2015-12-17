#!/usr/bin/env python

import os
import re
import math
import sys
import datetime
import time

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

delta = 10.0  # time of a fts spectrum , it migh be even longer

def match_fts_oco_times(oco_times_list, run_log_file, output_file):

    oco_times_fo = open(oco_times_list, 'r')

    oco_date_times = []
    oco_names = []
    for oco_time_str in oco_times_fo.readlines():
        oco_time_str = oco_time_str.strip()

        if not re.search('\d{16}', oco_time_str):
            continue
        
        oco_names.append(oco_time_str)

        oco_time_str = os.path.basename(oco_time_str).replace('.snd', '').replace('.dat', '')

        oco_year= long(oco_time_str[0:4])
        oco_mon  =long(oco_time_str[4:6])
        oco_day  =long(oco_time_str[6:8])
        oco_hour =long(oco_time_str[8:10])
        oco_min  =long(oco_time_str[10:12])
        oco_sec  =long(oco_time_str[12:14])
        oco_sidx =long(oco_time_str[14:15])
        #oco_frac = float(oco_h) + float(oco_m) / 60.0 + float(oco_s) / 3600.0 + long(oco_f) / 36000.0
        if oco_sidx == 0 or oco_sidx == 9:
            oco_fsec = 0
        elif oco_sidx == 3:
            oco_fsec = int(1e6/3.0)
        elif oco_sidx == 7 or oco_sidx == 6:
            oco_fsec = int(2e6/3.0)
        else:
            raise ValueError('Unknown frac seconding code: %d for time string: %s' % (oco_sidx, oco_time_str))
        
        oco_dt = datetime.datetime(oco_year, oco_mon, oco_day, oco_hour, oco_min, oco_sec, oco_fsec)
        oco_date_times.append(oco_dt)

    oco_times_fo.close()

    rl_fo = open(run_log_file, 'r')

    # Junk header line
    colnames = [ colname.lower() for colname in rl_fo.readline().split() ]

    fts_date_times = []
    fts_names = []
    for runlog_line in rl_fo.readlines():
        runlog_line = runlog_line.strip()
        runlog_parts = runlog_line.split()

        fts_name  = runlog_parts[colnames.index("spectrum_file_name")]
        fts_names.append(fts_name)
        
        fts_year  = long(runlog_parts[colnames.index("year")])
        fts_doy   = long(runlog_parts[colnames.index("day")])
        doy_dt    = datetime.datetime(*time.strptime('%d %d' % (fts_year, fts_doy), '%Y %j')[0:3])
        fts_hfrac = float(runlog_parts[colnames.index("hour")])
        
        if int(fts_hfrac) == 24:
            fts_hfrac = fts_hfrac - 24
            doy_dt += datetime.timedelta(1)
        
        fts_mon   = doy_dt.month
        fts_day   = doy_dt.day

        fts_hour = long(fts_hfrac)
        fts_minf = (fts_hfrac - fts_hour) * 60
        fts_min  = long(fts_minf)
        fts_fsec = (fts_minf - fts_min) * 60
        fts_sec  = long(fts_fsec)
        fts_fsec = long((fts_fsec - fts_sec) * 1e6)

        fts_dt = datetime.datetime(fts_year, fts_mon, fts_day, fts_hour, fts_min, fts_sec, fts_fsec)
        fts_date_times.append(fts_dt)
        
    rl_fo.close()

    num_fts = len(fts_names)
    out_obj = open(output_file, 'w')
    for (oco_name, oco_dt) in zip(oco_names, oco_date_times):
        for fts_idx in range(num_fts):
            fts_name = fts_names[fts_idx]
            fts_dt1  = fts_date_times[fts_idx]
            if fts_idx+1 < num_fts:
                fts_dt2 = fts_date_times[fts_idx+1]
            else:
                fts_dt2 = fts_dt1 + datetime.timedelta(seconds=delta)
                        
            if (oco_dt >= fts_dt1) and (oco_dt <= fts_dt2):
                print >>out_obj, oco_name, '\t', fts_name
                break
    out_obj.close()

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<oco_times_file> <run_log_file> <output_file>\n"
        sys.exit(1)

    oco_times_file  = sys.argv[1]
    run_log_file    = sys.argv[2]
    output_file     = sys.argv[3]

    match_fts_oco_times(oco_times_file, run_log_file, output_file)

if __name__ == "__main__":
    standalone_main()
