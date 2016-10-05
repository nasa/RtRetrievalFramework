#!/usr/bin/env python

from types import ListType
import re
import time
import os
import math
import sys

from OCO_MathUtil import *
from OCO_Matrix import OCO_Matrix

# this routine takes a runlog file and create a soundinginfo file
# this soundinginfo file can then be used for the convolution of FTS spectra
# to obtain OCO files with the correct header information

def create_sounding_info(runlog_file, sounding_info_file=None, spec_name_srch=None):

    runlog_fobj = open(runlog_file, "r")

    print 'Searching in runlog for:', spec_name_srch

    finished = False
    for runlog_line in runlog_fobj.readlines():
        runlog_parts = runlog_line.split()

        if (spec_name_srch == None or (spec_name_srch != None and re.search(spec_name_srch, runlog_parts[0]))):

            spec_name = runlog_parts[0]
            year = int(runlog_parts[1])
            doy  = int(runlog_parts[2])

            date_struct = time.strptime("%d %d" % (year, doy), "%Y %j")
            (month, day) = time.strftime("%m %d", date_struct).split()

            month = int(month)
            day   = int(day)
            
            frac     = float( runlog_parts[3] )
            hour     = long(frac)
            minute   = long((frac-float(hour))*60.0)
            sec      = long((frac-hour-float(minute)/60.0)*3600.0)
            sec_frac = (frac-hour-float(minute)/60.0)*3600.0-sec
  
            lat = runlog_parts[4]
            lon = runlog_parts[5]
            alt = runlog_parts[6]
            sza = runlog_parts[7]

            s_frac = '%5.4f' % sec_frac
            s_fr   = s_frac[1:5]

            timestamp=' frame_time_stamp = %04d-%02d-%02dT%02d:%02d:%02d%sZ' % (year, month, day, hour, minute, sec, s_fr)

            if sounding_info_file == None:
                sounding_info_file = 'soundinginfo_' + spec_name + '.dat'
            else:
                finished = True
                
            si_fobj = open(sounding_info_file, "w")

            si_fobj.write(timestamp + "\n")
            si_fobj.write(' sounding_altitude = ' + alt + "\n")
            si_fobj.write(' sounding_latitude = ' + lat + "\n")
            si_fobj.write(' sounding_longitude = ' + lon + "\n")
            si_fobj.write(' sounding_azimuth =  180.000000000000' + "\n")
            si_fobj.write(' sounding_zenith = ' + sza + "\n")
            si_fobj.close()
            
        if finished:
            break
            
    runlog_fobj.close()

def Process_File(fileObj, scriptOptions, valuesDict, mapDict):

    create_sounding_info(fileObj.filename, scriptOptions[0], scriptOptions[1])

def standalone_main():
    if (len(sys.argv) < 4):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<runlog_file> <sounding_input_filename> <spec_srch_string>\n"
        sys.exit(1)

    runlog_filename        = sys.argv[1]
    sounding_info_filename = sys.argv[2]
    spec_srch_string       = sys.argv[3]

    create_sounding_info(runlog_filename, sounding_info_filename, spec_srch_string)

if __name__ == "__main__":
    standalone_main()
