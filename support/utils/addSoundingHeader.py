#!/usr/bin/env python

from __future__ import print_function
import h5py
import sys
from numpy import *

import glob
from full_physics.l2_input import L2InputFile
import full_physics.acos_file as acos_file
from optparse import OptionParser

def addHeader(options):

    inp_config = L2InputFile(options.config_file)

    l1b_file = inp_config.get_section("input->InputProductFiles")[0].Get_Keyword_Value("L1BFile")
    print('Using '+  l1b_file + ' to extract header information')
    l1b_obj = acos_file.L1B(l1b_file)
    l2_obj = h5py.File(options.l2_file, "r+")
    try:
     #   print 'trying hard...'
        l2_obj.create_group("SoundingHeader")
        l2_obj.create_group("SoundingGeometry")
    except:
        print("Sounding header already exists in the L2 file!")
        sys.exit(0)
    exposure_index = l2_obj["RetrievalResults/exposure_index"][:]-1
    sounding_id_reference = l2_obj["/RetrievalResults/sounding_id_reference"][:]
    time_string = l1b_obj["SoundingHeader/exposure_start_time_string"][:][exposure_index]
    start_time = l1b_obj["SoundingHeader/exposure_start_time_tai93"][:][exposure_index]
    sounding_id = l1b_obj["SoundingHeader/sounding_id"][:][exposure_index]
    lat =  l1b_obj["SoundingGeometry/sounding_latitude"][:][exposure_index]
    lon = l1b_obj["SoundingGeometry/sounding_longitude"][:][exposure_index]
    sza = l1b_obj["SoundingGeometry/sounding_solar_zenith"][:][exposure_index]
    lza = l1b_obj["SoundingGeometry/sounding_zenith"][:][exposure_index]
    azi = l1b_obj["SoundingGeometry/sounding_solar_azimuth"][:][exposure_index]
    # write few entries into L2 files...
    l2_obj["SoundingHeader/exposure_start_time_string"] = time_string
    l2_obj["SoundingHeader/exposure_start_time_tai93"] = start_time
    l2_obj["SoundingHeader/sounding_id"] = sounding_id
    l2_obj["SoundingGeometry/sounding_longitude"] = lon
    l2_obj["SoundingGeometry/sounding_latitude"] = lat
    l2_obj["SoundingGeometry/sounding_solar_zenith"] = sza
    l2_obj["SoundingGeometry/sounding_zenith"] = lza
    l2_obj["SoundingGeometry/sounding_solar_azimuth"] = azi
    # sanity check:
    # print sounding_id_reference[0], sounding_id[0]
    l2_obj.close()
    l1b_obj.close()
    print('done...')

def standalone_main():
    parser = OptionParser(usage="usage: %prog --l2 l2_spliced.h5 --configFile configFile.config")

    parser.add_option( "--l2", dest="l2_file",
                       metavar="FILE",
                       help="spliced L2 dataset (or single sounding, whatever you like")
    parser.add_option( "--configFile", dest="config_file",
                       metavar="FILE",
                       help="config file used for the l2 run under consideration (the name of the L1b file is being pulled out of this file!)")
    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    if options.l2_file == None:
        parser.error('L2 file has to be specifiec')
    if options.config_file == None:
        parser.error('Config file has to be specified')
    addHeader(options)
        
    
if __name__ == "__main__":
    standalone_main()

#no guaranty, Christian Frankenberg
