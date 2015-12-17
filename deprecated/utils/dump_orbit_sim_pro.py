#!/usr/bin/env python

import os
import math
import sys
import struct

from OCO_MathUtil import *
from Orbit_Sim import Prof_File
from OCO_Matrix import OCO_Matrix

def dump_orbit_sim_pro(prof_file):
    prof_obj = Prof_File(prof_file)
     
    finished_file = False
    num_profiles = 0
    while not finished_file:
        if num_profiles != 0:
            print '=' * 75

        # Dump data
        profile_data = prof_obj.read_next_sounding(debug=True)
        if profile_data == None:
            finished_file = True

        num_profiles += 1
       
    prof_obj.close()
    
def Process_File(fileObj, scriptOptions, valuesDict, mapDict):
    pass

def standalone_main():
    if (len(sys.argv) < 2):
        print "usage:\n\t", os.path.basename(sys.argv[0]), "<orbit_sim_profile_file>\n"
        sys.exit(1)

    orbit_sim_profile_file  = sys.argv[1]

    dump_orbit_sim_pro(orbit_sim_profile_file)

if __name__ == "__main__":
    standalone_main()
