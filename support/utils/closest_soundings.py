#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import zip
from past.utils import old_div
import os
import sys
from optparse import OptionParser

import numpy
from pyproj import Geod

import full_physics.acos_file as acos_file

def find_closest_soundings(l1b_file, tgt_latitude, tgt_longitude, max_distance, log_output=sys.stdout):
    l1b_obj = acos_file.L1B(l1b_file)

    sounding_ids = l1b_obj.get_sounding_ids()

    # Average over band first since this is likely to be consistent for
    # different types of L1B files
    latitudes  = l1b_obj.get_sounding_info('sounding_latitude')
    longitudes = l1b_obj.get_sounding_info('sounding_longitude')

    # Average over any non sounding id sized dimensions
    while len(latitudes.shape) > 1:
        extra_dims = numpy.where(numpy.array(latitudes.shape) != sounding_ids.shape[0])
        latitudes  = numpy.average(latitudes, extra_dims[0][0])
        longitudes = numpy.average(longitudes, extra_dims[0][0])

    g = Geod(ellps='WGS84')

    # Find all distances in file
    distances = numpy.zeros(len(sounding_ids), dtype=float)
    for dist_idx, lat_lon_tuple in enumerate(zip(latitudes, longitudes)):
        curr_lat, curr_lon = lat_lon_tuple 
        az12, az21, dist = g.inv(tgt_longitude,tgt_latitude,curr_lon,curr_lat)
        
        # Convert to km
        distances[dist_idx] = old_div(dist,1000.)

    closest = numpy.where(distances <= max_distance)
    if len(closest[0]) > 0:
        print("%s" % l1b_file, file=log_output)
        for close_idx in closest[0]:
            print('%d %f' % (sounding_ids[close_idx], distances[close_idx]), file=log_output)
        print("", file=log_output)
    else:
        print("No soundings found in %s closer than %f km" % (l1b_file, max_distance), file=sys.stderr)
        
def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] --lat <latitude> --lon <longitude> <l1b_file>")

    parser.add_option( "--lat", dest="latitude",
                       type="float",
                       help="target latitude")

    parser.add_option( "--lon", dest="longitude",
                       type="float",
                       help="target longitude")

    parser.add_option( "-d", "--max_distance", dest="max_distance",
                       metavar="KM",
                       type="float",
                       default=100,
                       help="maximum distance in km between desired lat/lon and soundings, default is 100km")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('At least one l1b file needs to be specified')

    l1b_files = args

    if options.latitude == None:
        parser.error('Target latitude must be specified')

    if options.longitude == None:
        parser.error('Target longitude must be specified')

    print("Soundings less than %f km to latitude/longitude: %f, %f" % (options.max_distance, options.latitude, options.longitude), file=sys.stderr)

    for curr_l1b in l1b_files:
        find_closest_soundings(curr_l1b, options.latitude, options.longitude, options.max_distance)

if __name__ == "__main__":
    standalone_main()
