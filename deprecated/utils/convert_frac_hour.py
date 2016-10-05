#!/usr/bin/env python

import sys

frac = float( sys.argv[1] )
hours= long(frac)
mins = long((frac-float(hours))*60.)

sec = long((frac-hours-float(mins)/60.)*3600.)
sec_frac=(frac-hours-float(mins)/60.)*3600.-sec

print "%02d:%02d:%02d.%03d" % (hours, mins, sec, sec_frac)
