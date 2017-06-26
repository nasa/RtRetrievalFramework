#!/usr/bin/env python

from __future__ import print_function
from builtins import str
import h5py
import sys
import numpy as np
import glob
import datetime
import os
import stat
from optparse import OptionParser

# some global lists
orbit_list = ['004','008','012','018','019','023','027','033','037','041'];
#version_list = ['T_T_T','T_T_F','T_F_F','F_F_F'];
version_list = ['T_T_T','T_F_F','F_F_F'];

def setupRuns(options):
    # mainly intended to create the directory structure and set up the necessary l2 files...
    if not os.path.isdir(options.l2_directory):
        os.makedirs(options.l2_directory)
    # determine file name for postprocessing (shell script doing all necessary further steps)    
    postprocessingFile = options.l2_directory + '/postProcessing.sh'
    runFile = options.l2_directory + '/runL2.sh'
    if os.path.isfile(postprocessingFile):
        os.remove(postprocessingFile)
    if os.path.isfile(runFile):
        os.remove(runFile)    
    f = open(postprocessingFile,'w')
    r = open(runFile,'w')
        
    f.write('#!/bin/sh\n')
    f.write('# Shell script generated to automize postprocessing tools for the orbit simulator\n')
    f.write('# generated on ' + datetime.datetime.utcnow().strftime("%d/%m/%y, %HH:%MM") + '\n')
    f.write('# Christian Frankenberg\n')
    r.write('#!/bin/sh\n')
    r.write('# Shell script generated to automize running the l2 code\n')
    r.write('# generated on ' + datetime.datetime.utcnow().strftime("%d/%m/%y, %HH:%MM") + '\n')
    r.write('# Christian Frankenberg\n')

    # loop over orbits and run versions...
    for orbit in orbit_list:
        for version in version_list:
            # get name for the meteorology file, you may have to change foldernames to your needs here (most now set as input variable)!
            l1b = options.l1b_directory + options.l1b_subPath + '/7_' + version + '_orbit_' + orbit + '.hdf'
            met = options.l1b_directory + options.ecmwf_subPath +'/7_' + version +'_meteorology_' + orbit + '.hdf'
            logFile =  options.l1b_directory + options.l1b_subPath + '/7_' + version + '_orbit_' + orbit + '.log'
            sceneFile = options.l1b_directory + options.ecmwf_subPath +'/7_' + version +'_scene_' + orbit + '.log'
            hdf_logFile = options.l1b_directory + options.ecmwf_subPath + '/7_' + version +'_scene_' + orbit + '.hdf'
            opt_file = options.l1b_directory + options.l1b_subPath + '/7_' + version + '_orbit_' + orbit + '.opt'
            cloud_file = options.l1b_directory + options.cloud_subPath + version + '/7_' + version + '_Cld_' + orbit + '.hdf'
           
            
            conf = 'gosat_orbit_' + orbit + '_' + version + '.config'
            # check whether files exits and exit if not
            if not os.path.isfile(l1b) or not os.path.isfile(met):
                print('Matching l1b and/or met file not found: ' + l1b + met)
                break 
                #sys.exit(0)
           
            # determine and create the output directory
            outputDirectory = options.l2_directory + '/' + orbit + '/' + version
            if not os.path.isdir(outputDirectory):
                os.makedirs(outputDirectory)
            
            # run create_config.py to generate input file
            if not os.path.isfile(outputDirectory + '/' + conf) and not options.scriptOnly:
                print('Creating config file in ' + outputDirectory)
                if not options.ocean:
                    os.system('cd ' + outputDirectory + '; create_config.py -t gosat -g land -o ' + conf + ' ' + l1b + ' ' + met)
                else:
                    os.system('cd ' + outputDirectory + '; create_config.py -t gosat -g ocean -o ' + conf + ' ' + l1b + ' ' + met)
            else:
                if not options.scriptOnly:
                    print('Config file in ' + outputDirectory + ' already exists, doing nothing...')

            # run the populate tool to generate all necessary run directories, input files, etc.
            if not options.scriptOnly:
                print('Running the populate tool in ' + outputDirectory)
                os.system('cd ' + outputDirectory + '; populate.py -b ' + options.binary + '-q -o ' + conf)
        
            # for the time being: always rescale the co2 prior:
            if not options.scriptOnly:
                scale_CO2_prior(outputDirectory, 0.9873)
            if options.rayleighOnly:
                print("Scaling prior aerosol concentrations to ridicul. low values")
                os.system('sed -i \'s/-/-2/g\' input/static/scene/aerosol/profiles/aerosol_015_log_20.dat')
            
            # load constant ILS files
            if not options.scriptOnly:
                copyILS_files(os.path.abspath(outputDirectory))
                              
            # generate scripts for the l2 code runs
            r.write('# tools for running in ' +  outputDirectory + '\n')
            r.write('cd ' + os.path.abspath(outputDirectory) + '\n')
            r.write('./launch_jobs.sh -q amd \n')
            r.write('sleep 7000\n')

            # generate scripts for the postprocessing
            f.write('# tools for running in ' +  outputDirectory + '\n')
            f.write('cd ' + os.path.abspath(outputDirectory) + '\n')
            l2_final = 'l2_' + version + '_orbit_' +  orbit +  '.h5'
            f.write('splice_product_files.py -o '+ l2_final + ' output/*.h5\n')
            f.write('mergeOrbitSimOutput.py --l2='+l2_final+' --l1b='+l1b+' --logFile=' + logFile + ' --sceneFile=' + sceneFile + ' --detailedLog='+hdf_logFile+ ' --optFile='+opt_file+ ' --cloudFile='+cloud_file+'\n')
    r.close()
    f.close()
    
    os.chmod(postprocessingFile,stat.S_IXOTH)
    os.chmod(runFile,stat.S_IXOTH)
            
# ssh -o "LogLevel QUIET" fullerene 'cd orbitSim/benchmark; export HOSTNAME=`/bin/hostname` ; launch_jobs.py -r run_spec.dat -q long'            
#def postProcessing(options):

# Scales the CO2 prior profile
def scale_CO2_prior(dir, scaleFactor):
    print('Rescaling CO2 a priori profiles...')
    # find all matching CO2 a priori profiles:
    co2_files = glob.glob(dir+'/input/apriori/co2_apriori*.dat')
    for co2_file in co2_files:
        cmd = 'awk \'{if(NR>15) print $1,'+str(scaleFactor)+'*$2; else print $0}\' ' + co2_file + ' > ' + co2_file+'.bck'
        os.system(cmd)
        os.system('mv ' + co2_file+'.bck ' + co2_file) 
        #print co2_file

# Adapt this to your specific needs (where files are located...)
def copyILS_files(dir):
    print('Copying constant ILS files...')
    cmd = 'cp /home/cfranken/scratch/benchmark/7/ILS/band_1_v2.dat ' + dir+'/input/static/instrument/TFTSILSF_A006P006_B1S.dat'
    print(cmd)
    os.system(cmd)
    cmd = 'cp /home/cfranken/scratch/benchmark/7/ILS/band_2_v2.dat ' + dir+'/input/static/instrument/TFTSILSF_A006P006_B2.dat'
    print(cmd)
    os.system(cmd)
    cmd = 'cp /home/cfranken/scratch/benchmark/7/ILS/band_3_v2.dat ' + dir+'/input/static/instrument/TFTSILSF_A006P006_B3.dat'
    print(cmd)
    os.system(cmd)
   

# awk '{if(NR>15) print $1,$2; else print $0}' input/apriori/co2_apriori_20060914115552.dat

def standalone_main():
    parser = OptionParser(usage="usage: %prog \n This scripts performs several steps: \n 1) It creates level2 directories for all orbit simulator test cases \n 2) automatically runs create_config.py in all directories \n 3) automatically runs populate.py in all directories \n 4) copies constant ILS files into the respective rundirectories (as we ran benchmark 7 tests with constant ILS values) \n 5) Scales the CO2 prior profile with 0.973 (can be changed) \n 6) creates a shell script that automatically runs the l2 code in all subdirectories \n 7) creates a shell script that can automatically run the postprocessing mergin tool that generates single HDF5 files with all necessary information stored.")
    parser.add_option( "-d","--dir", dest="l1b_directory",
                       metavar="DIR",
                       default='/scratch2/validation/cfranken/benchmark/7/l1b_benchmark7_nadir/',
                       help="root directory where the l1b and met files are stored, default: /scratch2/validation/cfranken/benchmark/7/l1b_benchmark7_nadir/")
   
    parser.add_option( "-s","--scriptOnly", dest="scriptOnly",
                       action="store_true",default=False, help="only create pre and postprocessing shell scripts")
    parser.add_option( "-r","--rayleighOnly", dest="rayleighOnly",
                       action="store_true",default=False,help="set aerosol to zero in the L2 retrieval")
    parser.add_option( "--ocean", dest="ocean",
                       action="store_true",default=False,help="set this if yoy want to use only ocean pixels (otherwise only land pixels are used)")
    parser.add_option( "--l2_dir", dest="l2_directory",
                       metavar="DIR",
                       default='/scratch2/validation/cfranken/benchmark/7/l2_benchmark7_nadir',
                       help="root directory where the l2 results will be stored")
    parser.add_option( "--l1b_subPath", dest="l1b_subPath",
                       metavar="DIR",
                       default='/radiative_transfer/output/',
                       help="subpath to the l1b files from the l1b main directory, default:/radiative_transfer/output/")
    parser.add_option( "--ecmwf_subPath", dest="ecmwf_subPath",
                       metavar="DIR",
                       default='/scene_definition/output/',
                       help="subpath to the l1b files from the l1b main directory, default:/scene_definition/output/")
    parser.add_option( "--cloudScreen_subPath", dest="cloud_subPath",
                       metavar="DIR",
                       default='/cloudscreen/',
                       help="subpath to the l1b files from the l1b main directory, default:/cloudscreen/")
    parser.add_option( "-b","--binary", dest="binary", metavar='FILE',
                       default='/home/cfranken/acos/level_2/l2_fp',
                       help="Location of L2 binary to be used for runs")      
    # Parse command line arguments
    (options, args) = parser.parse_args()
    if options.scriptOnly:
        print("Just generating shell scripts for pre and postprocessing...") 
    setupRuns(options)
        
if __name__ == "__main__":
    standalone_main()

#no guarantee, Christian Frankenberg
