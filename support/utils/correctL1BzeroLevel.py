#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import sys
import h5py
from numpy import *
import numpy
from optparse import OptionParser

def correct_zero_level(curr_l1b, H_corr, M_corr, options):
    print(curr_l1b)
    f = h5py.File(curr_l1b,'r+')
    o2_spec = f["SoundingSpectra/radiance_o2"]
    gain = f["SoundingHeader/gain_swir"]
    try:
        f.create_group("ZeroLevelOffset")
        f["ZeroLevelOffset/GainMFile"] = options.m_file
        f["ZeroLevelOffset/GainHFile"] = options.h_file
        f["ZeroLevelOffset/GainMData"] = M_corr
        f["ZeroLevelOffset/GainHData"] = H_corr
    except:
        print("Datafile already corrected for zero level offset")
        sys.exit(0)
    
    for i in range(len(o2_spec[:,1,1])):
        sounding_gain = gain[i,0]
        average_signal_P =  o2_spec[i,0,:].mean()
        average_signal_S =  o2_spec[i,1,:].mean()
       
        if sounding_gain[0]=="H":
            corrP = 0.5*(numpy.interp(average_signal_P, H_corr[:,0], H_corr[:,1])+numpy.interp(average_signal_P, H_corr[:,0], H_corr[:,3]))
            corrS =0.5*(numpy.interp(average_signal_S, H_corr[:,0], H_corr[:,2])+numpy.interp(average_signal_S, H_corr[:,0], H_corr[:,4]))
            o2_spec[i,0,:] = o2_spec[i,0,:]-corrP*average_signal_P
            o2_spec[i,1,:] = o2_spec[i,1,:]-corrS*average_signal_S
         #   print average_signal_P, average_signal_S,  sounding_gain, corrP, corrS
        elif sounding_gain[0]=="M":
            corrP = 0.5*(numpy.interp(average_signal_P, M_corr[:,0], M_corr[:,1])+numpy.interp(average_signal_P, M_corr[:,0], M_corr[:,3]))
            corrS =0.5*(numpy.interp(average_signal_S, M_corr[:,0], M_corr[:,2])+numpy.interp(average_signal_S, M_corr[:,0], M_corr[:,4]))
            o2_spec[i,0,:] = o2_spec[i,0,:]-corrP*average_signal_P
            o2_spec[i,1,:] = o2_spec[i,1,:]-corrS*average_signal_S
          #  print average_signal_P, average_signal_S,  sounding_gain, corrP, corrS
    f.close()            


def standalone_main():
    parser = OptionParser(usage="usage: %prog -m M_gain_correction.dat -h H_gain_correction.dat l1B_files")

    parser.add_option( "--mFile", "--M_Gain_correctionFile", dest="m_file",
                       metavar="FILE",
                       help="text file with zero level offset corrections as a function of signal level (for Gain M)")
    parser.add_option( "--hFile", "--H_Gain_correctionFile", dest="h_file",
                       metavar="FILE",
                       help="text file with zero level offset corrections as a function of signal level (for Gain H)")
    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.error('At least one l1b file needs to be specified')
    if options.m_file == None:
        parser.error('M-Gain correction file must be specified')
    if options.h_file == None:
        parser.error('H-Gain correction file must be specified')
    l1b_files = args
    # load H and M gain correction curves...
    H_corr = numpy.loadtxt(options.h_file)
    M_corr = numpy.loadtxt(options.m_file)
    for curr_l1b in l1b_files:
        correct_zero_level(curr_l1b, H_corr, M_corr, options)
        
    
if __name__ == "__main__":
    standalone_main()
