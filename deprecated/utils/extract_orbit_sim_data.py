#!/usr/bin/env python

import os
import math
import sys
import struct
import logging
import traceback

from types import IntType
from optparse import OptionParser

import ACOS_File
from OCO_MathUtil import *
from Orbit_Sim import *
from OCO_Matrix import OCO_Matrix
import extract_acos_spectra
import L2_Log_Util

BAND_NAMES = ['abo2', 'wco2', 'sco2']
H2O_CONVERT_EPSILON = 0.621461

ALBEDO_CONVERSION = {
    'abo2': [0.33453167, -0.001975583],
    'wco2': [0.38253560, -0.003724288],
    'sco2': [0.53110406, -0.007071163],
    }

ALBEDO_CENTER_WAVELENGTHS = ( 0.77, 1.615, 2.06 )

LEV_BEG = 1

XCO2_COL_NAME     = 'XCO2 (ppmv)'
XCO2_LABEL_NAME   = 'TRUE_X_TARGET'
AOD_COL_NAMES = ['OT water', 'OT ice', 'OT aerosol']
AOD_LABEL_NAMES = ['WATER_OD', 'ICE_OD', 'AEROSOL_OD']
AOD_PROF_SRCH_LBL = ['wc_', 'ic_', ['cont', 'oceanic']]

DRY_AIR_NAME = 'AIR_dry'

def write_soundinginfo_file(hdf_sounding_dict, sounding_info_filename, sounding_id):
    sounding_info_fileobj = OCO_Matrix()
    sounding_info_fileobj.file_id = 'Sounding info from orbit simulator for sounding id: %s' % sounding_id

    sounding_info_fileobj.header = hdf_sounding_dict
    
    sounding_info_fileobj.header['sounding_id'] = sounding_id

    sounding_info_fileobj.write(sounding_info_filename)

def write_spectra_file(l1b_obj, hdf_sounding_dict, spectra_filename, sounding_id, average_name=None):
    spectra_fileobj = OCO_Matrix()
    extract_acos_spectra.write_ascii_from_hdf(l1b_obj, hdf_sounding_dict, spectra_fileobj, sounding_id, average_name=average_name)
    spectra_fileobj.write(spectra_filename, default_precision=16)
    
def write_xco2_file(log_sounding_dict, xco2_filename):
    xco2_fileobj = OCO_Matrix()
    xco2_fileobj.file_id = 'True xco2 from orbit simulator'
    xco2_fileobj.labels = [XCO2_LABEL_NAME]
    xco2_fileobj.data = numpy.zeros((1,1), dtype=float)
    xco2_fileobj.data[0,0] = log_sounding_dict[XCO2_COL_NAME]
    xco2_fileobj.write(xco2_filename)
    

def write_atmosphere_file(sounding_data, out_filename, pressure_in, pressure_out, profile_mass_densities=False):
    gas_names = sounding_data['gas_names']
    num_cols = len(gas_names) + 2
    
    out_atm_data = numpy.zeros((len(pressure_out), num_cols), dtype=float)

    out_atm_data[:, 0] = resample_profile(pressure_in, pressure_in, pressure_out, log_data=True, extrapolate=True)               
    out_atm_data[:, 1] = resample_profile(pressure_in, sounding_data['temperature'][LEV_BEG:], pressure_out, log_data=False, extrapolate=True)

    for gas_idx in range(len(gas_names)):
        curr_gas = sounding_data['gas_names'][gas_idx]

        # Ensure data is a numpy array 
        source_data = numpy.array(sounding_data['gas_mix_ratio'][curr_gas]) 

        if profile_mass_densities and curr_gas.find("AIR") != 0:
            # Convert H2O from mol/m^2 to VMR
            # The profiles in the Simulation structure are layer mass densities in mol/m^2
            # They are converted to volume mixing ratios by dividing by the mass density of dry air, also in mol/m^2
            source_data = source_data/sounding_data['gas_mix_ratio'][DRY_AIR_NAME]

        elif curr_gas == 'H2O':
            # Convert H2O from specific humidity%
            source_data = source_data/(1.0 - source_data)/H2O_CONVERT_EPSILON
                
        out_atm_data[:, gas_idx+2] = resample_profile(pressure_in, source_data[LEV_BEG:], pressure_out, log_data=True, extrapolate=False)
        

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'True atmospheric profiles from orbit simulator'
    out_mat_obj.data = out_atm_data
    out_mat_obj.labels = ['Pressure', 'T'] + list(gas_names)
    out_mat_obj.units =  ['Pa', 'K'] + [ 'VMR' for x in range(len(gas_names)) ]

    out_mat_obj.write(out_filename)

    return out_atm_data

def write_total_aod_file(log_sounding_dict, aod_filename):
    # make aerosol_od_<sounding_id>.dat
    
    aod_fileobj = OCO_Matrix()
    aod_fileobj.file_id = 'True aerosol optical depth from orbit simulator'
    aod_fileobj.labels = AOD_LABEL_NAMES
    aod_fileobj.data = numpy.zeros((1,len(AOD_COL_NAMES)), dtype=float)

    for out_idx, aer_col_name in enumerate(AOD_COL_NAMES):
        aod_fileobj.data[0,out_idx] = log_sounding_dict[aer_col_name]

    aod_fileobj.write(aod_filename)

def write_aerosol_file(sounding_data, out_filename, pressure_in, pressure_out):
    aer_names = sounding_data['aerosol_names']

    # If there are no aerosols for this sounding then just make an empty list
    if aer_names == None:
        aer_names = []
    
    num_cols = len(aer_names) + 1
   
    out_aer_data = numpy.zeros((len(pressure_out), num_cols), dtype=float)
   
    out_aer_data[:, 0] = resample_profile(pressure_in, pressure_in, pressure_out, log_data=True, extrapolate=True)               

    aer_src_data = numpy.zeros((len(pressure_in), len(aer_names)), dtype=float)
    
    for aer_idx in range(len(aer_names)):
        curr_aer = sounding_data['aerosol_names'][aer_idx]
        for lev_idx in range(len(pressure_in)):
            aer_src_data[lev_idx, aer_idx] = abs(sounding_data['aerosol_extinction'][curr_aer][LEV_BEG+lev_idx])

    if tuple(pressure_in) == tuple(pressure_out):
        aer_dst_data = aer_src_data.transpose()
    else:
        aer_dst_data = resample_aerosol(pressure_in, aer_src_data, pressure_out, debug=True)

    for aer_idx in range(len(aer_names)):
        out_aer_data[:, aer_idx+1] = aer_dst_data[aer_idx, :]

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'True aerosol profiles from orbit simulator'
    out_mat_obj.data = out_aer_data
    out_mat_obj.labels = ['Pressure'] + list(aer_names)
    out_mat_obj.header['Retrieval_Mode'] = 'linear'
    out_mat_obj.units =  ['Pa'] + [ '1/Pa' for x in range(len(aer_names)) ]

    out_mat_obj.write(out_filename)

def write_psurf_file(psurf, out_filename):

    out_psurf_data = numpy.zeros((1, 1), dtype=float)

    out_psurf_data[0, 0] = psurf

    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'True surface pressure from orbit simulator'
    out_mat_obj.data = out_psurf_data
    out_mat_obj.labels = ['PSURF']
    out_mat_obj.units =  ['Pa']

    out_mat_obj.write(out_filename)

def write_windspeed_file(sounding_data, out_filename):

    out_windspeed_data = numpy.zeros((1, 1), dtype=float)

    ws_data = sounding_data['surface_windspeed']
    if hasattr(ws_data, '__iter__'):
        out_windspeed_data[0, 0] = ws_data[0]
    else:
        out_windspeed_data[0, 0] = ws_data
        
    out_mat_obj = OCO_Matrix()
    out_mat_obj.file_id = 'True windspeed from orbit simulator'
    out_mat_obj.data = out_windspeed_data
    out_mat_obj.labels = ['WINDSPEED']
    out_mat_obj.units =  ['m/s']

    out_mat_obj.write(out_filename)

def write_albedo_file(log_sounding_dict, out_filename):

    out_albedo_data = numpy.zeros((2, len(BAND_NAMES)), dtype=float)

    for band_idx in range(len(BAND_NAMES)):
        albedo_os_a = log_sounding_dict['EffAlbedo_%dA' % (band_idx+1)]
        albedo_os_b = log_sounding_dict['EffAlbedo_%dB' % (band_idx+1)]

        conv_f = ALBEDO_CONVERSION[BAND_NAMES[band_idx]][0]
        conv_g = ALBEDO_CONVERSION[BAND_NAMES[band_idx]][1]
        
        albedo_l2 = albedo_os_a * (1-conv_f) + albedo_os_b * conv_f
        slope = (albedo_os_b - albedo_os_a) * conv_g

        out_albedo_data[0, band_idx] = albedo_l2
        out_albedo_data[1, band_idx] = slope

    out_mat_obj = OCO_Matrix()
    out_mat_obj.header['center_wavelengths'] = ' '.join([str(wl) for wl in ALBEDO_CENTER_WAVELENGTHS])
    out_mat_obj.file_id = 'True albedo from orbit simulator'
    out_mat_obj.data = out_albedo_data
    out_mat_obj.labels = ['ALBEDO_%d' % (bidx+1) for bidx in range(out_albedo_data.shape[1])]

    out_mat_obj.write(out_filename)


def extract_orbit_sim_data(prof_file, l1b_file, true_dir, log_files=None, resample_to=None, sounding_id_list=None, group_sounding_data=False, filter_retrievable=False, write_spectra=None, verbose=False):
    logger = logging.getLogger(os.path.basename(__file__))

    if not os.path.exists(true_dir):
        raise IOError('true destination directory does not exist: %s' % true_dir)

    # Load sounding list
    if not os.path.exists(l1b_file):
        raise IOError('L1B file specified does not exist: %s' % l1b_file)
    
    l1b_obj = ACOS_File.L1B(l1b_file)
    snd_id_matrix = l1b_obj.get_sounding_ids()
   
    # Open input file in binary mode
    if not os.path.exists(prof_file):
        raise IOError('profile file does not exist: %s' % prof_file)

    if log_files != None and len(log_files) != 0:
        log_file_objs = [ Log_File(curr_file) for curr_file in log_files ]
    else:
        log_file_objs = None

    profile_mass_densities = False
    if h5py.is_hdf5(prof_file):
        # Newer GOSAT simulator profile information file
        prof_obj = Simulation_File(prof_file)
        profile_mass_densities = True
    else:
        # Older OCO simulator profile information file
        prof_obj = Prof_File(prof_file)
   
    if filter_retrievable and log_file_objs != None:
        performed_filter = False
        for curr_log_obj in log_file_objs:
            try:
                if sounding_id_list != None:
                    sounding_id_list = curr_log_obj.filter_retrievable(sounding_id_list)
                else:
                    sounding_id_list = curr_log_obj.filter_retrievable(numpy.ravel(snd_id_matrix))
                performed_filter = True
            except LookupError:
                print >>sys.stderr, 'Could not use %s for filtering sounding' % curr_log_obj.filename
                for error_line in traceback.format_exception(*sys.exc_info()):
                    print >>sys.stderr, error_line


        if not performed_filter:
            raise IOError('Could not filter sounding list using any available log file: %s' % ([log_obj.filename for log_obj in log_file_objs]))
                

    if resample_to != None:
        if resample_to.isdigit():
            resample_to = int(resample_to)
        else:
            pres_obj = OCO_Matrix(resample_to)
            src_pres_col = pres_obj.labels_lower.index("pressure")
            resample_to = pres_obj.data[:, src_pres_col]

    if l1b_obj.instrument_name == ACOS_File.GOSAT_INST_NAME:
        average_name = 'Polarization'
    else:
        average_name = None

    for pro_sounding_data in prof_obj:
        # Determine sounding id based on current index in prof_obj
        if isinstance(prof_obj, Prof_File):
            frame_index    = pro_sounding_data['frame_index']
            sounding_index = pro_sounding_data['sounding']-1

            sounding_id = long(snd_id_matrix[frame_index, sounding_index])
        elif isinstance(prof_obj, Simulation_File):
            sounding_id = long(snd_id_matrix[ pro_sounding_data.exposure_index ])            
        else:
            raise Exception('prof_obj is an unknown object: %s' % prof_obj)

        # Make sure the current sounding is in the list of those to process
        if sounding_id_list != None and not sounding_id in sounding_id_list:
            continue

        if verbose:
            logger.info('Processing sounding: %s' % sounding_id)

        if group_sounding_data:
            output_dir = os.path.join(true_dir, str(sounding_id))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        else:
            output_dir = true_dir

        # Get dictionary of values from all log files 
        if log_file_objs != None:
            log_sounding_dict = {}
            for curr_log_obj in log_file_objs:
                log_sounding_dict = dict(curr_log_obj.get_sounding_info_dict(sounding_id), **log_sounding_dict)

        # Load from the profile object the pressure column for use in interpolation when necessary
        pres_in = pro_sounding_data['pressure'][LEV_BEG:]
        if resample_to == None:
            pres_out = pres_in
        else:
            pres_out = resample_profile(pres_in, pres_in, resample_to, log_data=True, extrapolate=True)


        # Get sounding dictionary for ascii spectra header and sounding info file
        hdf_sounding_dict = l1b_obj.get_sounding_info_dict(sounding_id, ignore_missing=True, average=average_name)
        soundinginfo_true_filename = os.path.join(output_dir, 'soundinginfo_%s.dat' % sounding_id)
        if verbose:
            logger.info('Writing sounding info file: %s' % soundinginfo_true_filename)
            write_soundinginfo_file(hdf_sounding_dict, soundinginfo_true_filename, sounding_id)
        
        # If L1B is a GOSAT file then also write the per polarization information as well as averaged information from above
        # Create sounding file from L1B information using a dictionary of values from L1B file
        sounding_head_ids = []
        if l1b_obj.instrument_name == ACOS_File.GOSAT_INST_NAME:
            for pol_name in ACOS_File.GOSAT_POL_ORDER:
                sounding_id_pol = '%s%s' % (sounding_id, pol_name)
                sounding_head_ids.append(sounding_id_pol)

        for curr_head_id in sounding_head_ids:
            pol_sounding_dict = l1b_obj.get_sounding_info_dict(sounding_id, ignore_missing=True)
        
            soundinginfo_true_filename = os.path.join(output_dir, 'soundinginfo_%s.dat' % curr_head_id)
            if verbose:
                logger.info('Writing sounding info file: %s' % soundinginfo_true_filename)
            write_soundinginfo_file(pol_sounding_dict, soundinginfo_true_filename, curr_head_id)

        if log_file_objs != None: 
            xco2_true_filename = '%s/xco2_%d.dat' % (output_dir, sounding_id)
            if verbose:
                logger.info('Writing xco2 file: %s' % xco2_true_filename)
            write_xco2_file(log_sounding_dict, xco2_true_filename)

        atm_true_filename = '%s/atmosphere_%d.dat' % (output_dir, sounding_id)
        if verbose:
            logger.info('Writing atmosphere file: %s' % atm_true_filename)
        write_atmosphere_file(pro_sounding_data, atm_true_filename, pres_in, pres_out, profile_mass_densities)

        if log_file_objs != None:
            aod_true_filename = '%s/aerosol_od_%d.dat' % (output_dir, sounding_id)
            if verbose:
                logger.info('Writing aerosol optical depth file: %s' % aod_true_filename)
            write_total_aod_file(log_sounding_dict, aod_true_filename)

        if l1b_obj.instrument_name == ACOS_File.OCO_INST_NAME:
            # Only write these for OCO since the translation of units
            # for GOSAT simulator files DOES NOT yet work
            aer_true_filename = '%s/aerosol_%d.dat' % (output_dir, sounding_id)
            if verbose:
                logger.info('Writing aerosol file: %s' % aer_true_filename)
            write_aerosol_file(pro_sounding_data, aer_true_filename, pres_in, pres_out)

        psurf_true_filename = '%s/psurf_%d.dat' % (output_dir, sounding_id)
        if verbose:
            logger.info('Writing psurf file: %s' % psurf_true_filename)
        write_psurf_file(pres_in[-1], psurf_true_filename)

        windspeed_true_filename = '%s/windspeed_%d.dat' % (output_dir, sounding_id)
        if verbose:
            logger.info('Writing windspeed file: %s' % windspeed_true_filename)
        write_windspeed_file(pro_sounding_data, windspeed_true_filename)

        if log_file_objs != None:
            albedo_true_filename = '%s/albedo_%d.dat' % (output_dir, sounding_id)
            if verbose:
                logger.info('Writing albedo file: %s' % albedo_true_filename)
            write_albedo_file(log_sounding_dict, albedo_true_filename)

        if write_spectra != None:
            spectra_filename = '%s/spectra_%d.dat' % (output_dir, sounding_id)
            if verbose:
                logger.info('Writing spectra file: %s' % spectra_filename)
            write_spectra_file(l1b_obj, hdf_sounding_dict, spectra_filename, sounding_id, average_name)

    prof_obj.close()
       
if __name__ == "__main__":    
    parser = OptionParser(usage="usage: %prog [options] <prof_file> <l1b_file> [<log_file>]\n")
    
    parser.add_option( "-t", "--true_dir", dest="true_dir",
                       metavar="DIR", default="./",
                       help="location to output true files",
                       )

    parser.add_option( "-g", "--group_sounding", dest="group_sounding_data",
                       action="store_true",
                       help="make a subdir under true_dir named for sounding id",
                       )

    parser.add_option( "-f", "--filter_retrievable", dest="filter_retrievable",
                       action="store_true",
                       help="filter out only retievable soundings for extraction",
                       )

    parser.add_option( "-l", "--log_file", dest="log_files",
                       metavar="FILE",
                       action='append',
                       help="location orbit simulator log file",
                       )

    parser.add_option( "-r", "--resample_to", dest="resample_to",
                       metavar="INT_OR_FILE",
                       help="destination pressure grid file or number of destination levels",
                       )

    parser.add_option( "-s", "--sounding_list", dest="sounding_list_file",
                       metavar='FILE',
                       help="sounding id list file specifiying which to extract"
                       )

    parser.add_option( "--spectra", dest="write_spectra",
                       action="store_true",
                       help="write spectra files from hdf file",
                       )


    parser.add_option( "-v", "--verbose", dest="verbose",
                       action="store_true",
                       help="enable verbose output information",
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if (len(args) < 2):
        parser.error("All required arguments not specified")
    
    prof_file = args[0]
    l1b_file  = args[1]

    L2_Log_Util.init_logging()

    sounding_ids = None
    if options.sounding_list_file != None:
        with open(options.sounding_list_file, 'r') as snd_list_obj:
            sounding_ids = [ long(line.strip()) for line in snd_list_obj.readlines() ]

    extract_orbit_sim_data(prof_file, l1b_file, log_files=options.log_files, true_dir=options.true_dir, resample_to=options.resample_to, group_sounding_data=options.group_sounding_data, filter_retrievable=options.filter_retrievable, write_spectra=options.write_spectra, verbose=options.verbose, sounding_id_list=sounding_ids)
