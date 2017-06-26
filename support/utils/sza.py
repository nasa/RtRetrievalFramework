#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import datetime
import time
import math
import sys
import os
import re

from optparse import OptionParser

# PyEphem, available at http://www.rhodesmill.org/brandon/projects/pyephem.html
import ephem

def computeSZA(sun, location):
    sun.compute(location)
    return (old_div(math.pi,2) - sun.alt)

def bisect(sun, location, start, end, sza):
    x0 = start
    x1 = end

    orig_date = location.date

    location.date = x0
    y0 = computeSZA(sun, location) - sza

    location.date = x1
    y1 = computeSZA(sun, location) - sza
    
    if (y0 > y1):
        tmp = x0
        x0 = x1
        x1 = tmp

    while (math.fabs(x1 - x0) > 1e-6):
        mid = old_div((x0+x1), 2)
        location.date = mid

        curr_sza = computeSZA(sun, location)

        if (curr_sza - sza < 0):
            x0 = mid
        else:
            x1 = mid

    return mid

def getTimeOfSZA(lat, lon, date, desired_sza=None):
    r2d = old_div(180,math.pi)
    tolerance = 0.5 # if min or max sza is within this range of sza, accept it

    success = True

    location = ephem.Observer()
    location.lat = lat
    location.long = lon
    location.date = date

    sun = ephem.Sun()
    sun.compute(location)

    if desired_sza != None:
        noon = sun.transit_time

        if (noon == None):
            print("No transit at latitude %s, longitude %s at %s" % (lat,lon,date))
            return (location.date, 999, False)

        midnight = noon + ephem.hour * 12
        location.date = midnight
        max_sza = r2d*computeSZA(sun, location)

        if (max_sza - desired_sza > -tolerance):
            location.date = noon
            min_sza = r2d*computeSZA(sun, location)

            if (desired_sza - min_sza > -tolerance):
                location.date = bisect(sun, location, date - ephem.hour * 12, date + ephem.hour * 12, old_div(desired_sza,r2d))
                success = True
            else:
                # The sun never gets high enough to reach desired sza
                print("Minimum sza at %s is %f, desired is %f, difference is %f" % (dateToString(location.date), min_sza, desired_sza, desired_sza - min_sza))
                success = False
        else:
            # The sun never gets low enough to reach sza (polar day)
            print("Maximum sza at %s is %f, desired is %f, difference is %f" % (dateToString(location.date), max_sza, desired_sza, max_sza - desired_sza))
            success = False

    return (location.date, r2d*computeSZA(sun, location), success)

def readFile(filename, keywords):
    returnval = [None for x in keywords]
    try:
        file = open(filename, "r")

        while(True):
            this_line = file.readline()
            if (len(this_line) == 0): break

            # Skip commented lines
            hash = this_line.find("#")
            if (hash != -1):
                this_line = this_line[:hash]
            if (len(this_line) == 0): continue
                
            line = this_line.lower()

            tokens = line.split()
            if (len(tokens) == 0): continue

            for i in range(len(keywords)):
                if (tokens[0] == keywords[i]):
                    returnval[i] = tokens[2]
            
        file.close()
    except IOError:
        print("%s does not exist" % filename)
        sys.exit()

    return returnval

# Replace all entries with specified keywords with corresponding values
def replaceEntry(filename, keywords, values):
    try:
        file = open(filename, "r")
        inp_lines = file.readlines()
        file.close()
            
        newfile = open("%s" % filename, "w")

        for this_line in inp_lines:
            line = this_line.lower()
            tokens = line.split()
            if (len(tokens) == 0): continue

            for i in range(len(keywords)):
                if (tokens[0] == keywords[i]):
                    orig_tokens = this_line.split()
                    this_line = "%s = %s\n" % (orig_tokens[0], values[i])
            newfile.write(this_line)
            
        file.close()
        newfile.close()
    except IOError:
        print("%s does not exist" % filename)
        sys.exit()

    return

# take a ephem.date and convert to a UTC string
def dateToString(date):
    string = "%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%06.3fZ" % date.tuple()
    return string

def dateToFTS(date):
    (year, month, day, hour, minute, second) = date.tuple()

    doy = time.strftime('%j', datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second)).timetuple())

    frac_hour = hour + (old_div(minute, 60.0)) + (old_div(second, 3600.0))

    string = "%4.4d %03.3d %6.3f" % (int(year), int(doy), float(frac_hour))
    return string

# take a UTC string and create an ephem date
def stringToDate(UTC):
    year = int(UTC[0:4])
    month = int(UTC[5:7])
    day = int(UTC[8:10])
    hour = int(UTC[11:13])
    minute = int(UTC[14:16])
    second = float(UTC[17:23])

    date = ephem.date((year, month, day, hour, minute, second))
    return date

def getUsage():
    usage = ""
    usage += "usage:\n\t%s soundinfo_file [options]\n" % os.path.basename(sys.argv[0])
    usage += "e.g. %s sounding_info.dat\n" % os.path.basename(sys.argv[0])
    usage += "or\n\t%s lat lon date [options]\n" % os.path.basename(sys.argv[0])
    usage += "e.g. %s 46 -90.25 2007-03-20T00:00:00.000Z\n" % os.path.basename(sys.argv[0])
    
    return usage

if __name__ == "__main__":

    # Load command line options
    parser = OptionParser(usage=getUsage())
    
    parser.add_option( "-s", "--desired_sza", dest="desired_sza",
                       default=None, 
                       help="desired solar zenith angle around given time" )

    parser.add_option( "-f", "--fts_time", dest="fts_time",
                       default=False,
                       action="store_true",
                       help="print using fts time" )

    # Parse command line arguments
    arg_process = []
    args = []
    p_next = False
    for curr_a in sys.argv[1:]:
        if re.match('-[a-zA-Z].*', curr_a):
            arg_process.append(curr_a)
            p_next = True
        elif p_next:
            arg_process.append(curr_a)
            p_next = False
        else:
            args.append(curr_a)
            
    (options, parsed_args) = parser.parse_args(args=arg_process)

    
    if options.desired_sza != None:
        options.desired_sza = float(options.desired_sza)

    if (len(args) == 1):
        sounding_file = args[0]

        keywords = ["frame_time_stamp", "sounding_latitude", "sounding_longitude"]
        returnval = readFile(sounding_file, keywords)
        lat = returnval[1]
        lon = returnval[2]
        date = stringToDate(returnval[0])

    elif (len(args) == 3):
        lat = args[0]
        lon = args[1]
        date = stringToDate(args[2])
        
    else:
        parser.error('unrecognized number of arguments')

        
    returnval = getTimeOfSZA(lat, lon, date, options.desired_sza)

    if returnval[2]:
        if options.fts_time:
            print("sza at %s is %f" % (dateToFTS(returnval[0]), returnval[1]))
        else:
            print("sza at %s is %f" % (dateToString(returnval[0]), returnval[1]))
    else:
        raise ValueError("error calculating sza")
