#!/usr/bin/env python

# Load standard modules
import os
import sys
import copy
import glob
import bisect
import logging
import traceback
from optparse import OptionParser

# Load user installed modules
import numpy

# Load L2 modules
from OCO_TextUtils import *
import OCO_MathUtil
from OCO_Matrix import OCO_Matrix
import L2_Log_Util

##
LOGGING_NAME = 'testset_summary'

## For summary file

FILL_VALUE = None

# Used for AOD calculation
ATMOSPHERE_TRUE_BASE = 'true/atmosphere*.dat'
PRESSURE_TRUE_BASE = ATMOSPHERE_TRUE_BASE
PRESSURE_RET_BASE  = 'out/atmosphere.dat.spec??'

NEW_LEVEL_FRAC     = .20

PRESSURE_COLUMN    = 'Pressure'
TEMPERATURE_COLUMN = 'T'
H2O_COLUMN         = 'H2O'
CO2_COLUMN         = 'CO2'
O2_COLUMN          = 'O2'

PSURF_TRUE_BASE    = 'true/psurf*.dat'
PSURF_COLUMN       = 'PSURF'

VALID_PSURF_MIN    = 400.00

RESULT_FILE_BASE = 'out/results.dat'
STATEVECTOR_BASE = 'out/control1/final_statevector.dat'

RAD_MEAS_BASE    = 'out/rad_meas.dat'
RADIANCE_COLUMN  = 'Radiance'
ERROR_COLUMN     = 'Error'
WAVELENGTH_COLUMN = 'Wavelength'

SV_NAMES_BASE    = 'out/aggregator/sv_names.dat'
SV_NAMES_COLUMN  = 'Element Name'

LARGEST_PRESSURE = 1.0e20

#  Pressure ranges evaluted:     > (gt)   <= (le)
AOD_SUB_TYPES = { 'TOTAL_AOD': (0,     LARGEST_PRESSURE),
                  'LOW_AOD'  : (800e2, LARGEST_PRESSURE),
                  'MID_AOD'  : (500e2, 800e2),
                  'HIGH_AOD' : (0,     500e2),
                  }
AOD_TYPE_ORDER = ( 'TOTAL_AOD', 'LOW_AOD', 'MID_AOD', 'HIGH_AOD' )
AOD_CHECK_EVAL = 'abs({TOTAL_AOD} - ({LOW_AOD} + {MID_AOD} + {HIGH_AOD})) < 1e-6'

PER_RUN_CHECKS = ( ( ('CALC_TRUE_ALL_TOTAL_AOD', 'CALC_TRUE_ICE_TOTAL_AOD', 'CALC_TRUE_WATER_TOTAL_AOD', 'CALC_TRUE_CLAT_0_TOTAL_AOD', 'CALC_TRUE_CLAT_1_TOTAL_AOD', 'CALC_TRUE_CLAT_2_TOTAL_AOD', 'CALC_TRUE_CLAT_3_TOTAL_AOD', 'CALC_TRUE_CLAT_4_TOTAL_AOD', 'CALC_TRUE_CLAT_5_TOTAL_AOD', 'CALC_TRUE_CLAT_6_TOTAL_AOD'),
                     'abs({CALC_TRUE_ALL_TOTAL_AOD} - ({CALC_TRUE_ICE_TOTAL_AOD} + {CALC_TRUE_WATER_TOTAL_AOD} + {CALC_TRUE_AERO1_TOTAL_AOD} + {CALC_TRUE_AERO2_TOTAL_AOD} + {CALC_TRUE_CLAT_0_TOTAL_AOD} + {CALC_TRUE_CLAT_1_TOTAL_AOD} + {CALC_TRUE_CLAT_2_TOTAL_AOD} + {CALC_TRUE_CLAT_3_TOTAL_AOD} + {CALC_TRUE_CLAT_4_TOTAL_AOD} + {CALC_TRUE_CLAT_5_TOTAL_AOD} + {CALC_TRUE_CLAT_6_TOTAL_AOD})) < 1e-6'),
                     
                   ( ('RET_ALL_TOTAL_AOD', 'RET_ICE_TOTAL_AOD', 'RET_WATER_TOTAL_AOD', 'RET_CONT_TOTAL_AOD', 'RET_OCEANIC_TOTAL_AOD'),
                     'abs({RET_ALL_TOTAL_AOD} - ({RET_ICE_TOTAL_AOD} + {RET_WATER_TOTAL_AOD} + {RET_CONT_TOTAL_AOD} + {RET_OCEANIC_TOTAL_AOD})) < 1e-9' ),
                   )

RUN_COLUMN_NAME = 'Case'

SV_RESULT_COL  = 'Statevector'
SV_APRIORI_COL = 'A_Priori'

# Corrections for problematic strings when running extract result values
EXTRACT_RESULT_CORRECTIONS = ( ('-999$', ' -999'), )  # Fix missing space before -999 in some results.dat files

# Extract routines placed here since need to be known before CONFIG dicts
def extract_file_value(ts_obj, file_obj, curr_dir, row_idx, column_desc, data_label):

    if type(column_desc) == int:
        column_indexes = [ column_desc ]
    else:
        column_indexes = file_obj.find_labels(column_desc, indexes=True)

    if len(column_indexes) == 0:
        raise LookupError('Could not find column: %s in file %s' % (column_desc, file_obj.filename))
    elif len(column_indexes)  > 1:
        raise LookupError('Too many columns found named: %s in file %s'  % (column_desc, file_obj.filename))

    try:
        file_value = file_obj[column_indexes[0]][row_idx]
    except IndexError:
        raise IndexError('Invalid index %d, %d in file %s with shape %s' % (row_idx, column_indexes[0], file_obj.filename, file_obj.data.shape))
    
    return [ ( data_label, file_value ) ]

def calculate_file_xco2(ts_obj, file_obj, curr_dir, psurf_file, ak_file, data_label):
    logger = logging.getLogger(LOGGING_NAME)

    # See if extra files are available
    psurf_obj = ts_obj.get_cached_file(curr_dir, psurf_file)
    ak_obj = ts_obj.get_cached_file(curr_dir, ak_file)

    if ak_obj == None:
        ak_data = None
    else:
        ak_data = ak_obj.data

    try:
        in_press_profile = file_obj[PRESSURE_COLUMN]
    except LookupError:
        logger.error('For file "%s" LookupError: %s' % (file_obj.filename,traceback.format_exception_only(*sys.exc_info()[0:2])))
        return None

    if in_press_profile == None or in_press_profile.shape[0] == 0:
        logger.error('Could not calculate xco2 from file: "%s" for run dir "%s" due to missing pressure column' % (file_obj.filename. curr_dir))

    if psurf_obj != None:
        out_press_profile = apply_psurf_to_pressures(in_press_profile[:,0], psurf_obj[PSURF_COLUMN][0,0])
    else:
        out_press_profile = in_press_profile[:,0]


    # Get desired atmosphere profiles for calculation
    try:
        co2_profile = file_obj[CO2_COLUMN]
    except LookupError:
        logger.error('Could not find column named: %s' % CO2_COLUMN)
        return None
    
    if co2_profile == None or co2_profile.shape[0] == 0:
        logger.error('Could not calculate xco2 from file: "%s" for run dir "%s" due to missing CO2 column' % (file_obj.filename. curr_dir))
        return None
    
    try:
        co2_profile = OCO_MathUtil.resample_profile(in_press_profile, co2_profile[:,0], out_press_profile)
    except ValueError:
        logger.error('Could not resample profile with inputs: in_pressure_profile = %s, co2_profile = %s, out_pressure_profile = %s' % (in_press_profile, co2_profile[:,0], out_press_profile))
        return None

    # Get temperature and water profiles, if results are empty make sure we
    # pass None to compute_xco2
    temp_profile = file_obj[TEMPERATURE_COLUMN]
    if temp_profile != None and temp_profile.shape[0] == 0:
        temp_profile = None
    else:
        temp_profile = OCO_MathUtil.resample_profile(in_press_profile, temp_profile[:,0], out_press_profile)

    h2o_profile = file_obj[H2O_COLUMN]
    if h2o_profile != None and h2o_profile.shape[0] == 0:
        h2o_profile = None
    else:
        h2o_profile = OCO_MathUtil.resample_profile(in_press_profile, h2o_profile[:,0], out_press_profile)

    calc_xco2 = OCO_MathUtil.compute_xco2(co2_profile, out_press_profile, temp_profile, h2o_profile, ak_data)

    return [ (data_label, calc_xco2) ]

def extract_header_value(ts_obj, file_obj, curr_dir, keyword_name, data_label):
    logger = logging.getLogger(LOGGING_NAME)

    if not file_obj.header.has_key(keyword_name):
        logger.error('Header of "%s" does not contain keyword "%s"' % (file_obj.filename, keyword_name))
        return None

    return [ (data_label, file_obj.header[keyword_name]) ]

def extract_result_value(ts_obj, file_obj, curr_dir, row_idx, column_desc, data_label):
    logger = logging.getLogger(LOGGING_NAME)
    
    # Reorder structure of result file as read from OCO_Matrix
    if not hasattr(file_obj, 'fixed_results'):
        if not type(file_obj.data) is list:
            raise Exception('Expected result file: %s to be parsed as a list of lists not: %s' % (file_obj.filename, type(file_obj.data)))
        
        new_lines = []
        curr_line = []
        num_new_cols = 0
        for row_idx in range(len(file_obj.data)):
            in_row_data = file_obj.data[row_idx]

            # Fix row data
            fixed_row_data = []
            for row_value in in_row_data:
                for from_str, to_str in EXTRACT_RESULT_CORRECTIONS:
                    row_value = re.sub(from_str, to_str, row_value)
                    fixed_row_data += row_value.split()

            if fixed_row_data[0].find('Result') >= 0:
                if len(curr_line) > 0:
                    new_lines.append(curr_line)
                    num_new_cols = max(num_new_cols, len(curr_line))
                    curr_line = []
                num_new_cols = max(num_new_cols, len(fixed_row_data))
                new_lines.append(list(fixed_row_data))
            else:
                curr_line += fixed_row_data
                
        if len(curr_line) > 0:
            new_lines.append(curr_line)
            num_new_cols = max(num_new_cols, len(curr_line))
            curr_line = []


        new_matrix = numpy.zeros((len(new_lines), num_new_cols), dtype=numpy.chararray)
        for row_idx, row_data in enumerate(new_lines):
            for col_idx, col_val in enumerate(row_data):
                new_matrix[row_idx, col_idx] = col_val

        file_obj.data = new_matrix
        file_obj.dims = new_matrix.shape
        file_obj.fixed_results = True

    if file_obj.dims[0] <= 1:
        logger.error('Invalid results file: %s for run directory: %s' % (file_obj.filename, curr_dir))
        return None

    return extract_file_value(ts_obj, file_obj, curr_dir, row_idx, column_desc, data_label)

def extract_outinfo_value(ts_obj, file_obj, curr_dir, column_desc, data_label):
    logger = logging.getLogger(LOGGING_NAME)
    
    # Reorder structure of result file as read from OCO_Matrix
    if not hasattr(file_obj, 'fixed_outinfo'):
        if not type(file_obj.data) is list:
            raise Exception('Expected out_info file: %s to be parsed as a list of lists not: %s' % (file_obj.filename, type(file_obj.data)))

        new_lines = []
        file_labels = []
        file_obj.data.reverse()
        for row_data in file_obj.data:
            if row_data[0].find('kat') >= 0 and len(new_lines) > 0:
                nspec = int( (len(new_lines[-1]) - 10) / 2 )
                for curr_label in row_data:
                    if curr_label.find('(1:nspec)') >= 0:
                        curr_label = curr_label.replace('(1:nspec)', '')
                        for idx in range(nspec):
                            file_labels.append( '%s_%d' % (curr_label, idx+1) )
                    elif curr_label.find('error_flag') >= 0 or curr_label == 'e':
                        for idx in range(2):
                            file_labels.append( '%s_%d' % ('error_flag', idx+1) )
                    else:
                        file_labels.append(curr_label)
                file_obj.labels = ['orbit'] + file_labels
                break
            elif not re.search(row_data[0], '\d+') and len(new_lines) > 0:
                # Ignore wrapped column that has no "data" aka numbers
                pass
            else:
                new_lines.insert(0, row_data)

        file_obj.fixed_outinfo = True

        if len(new_lines) > 1:
            file_obj.data = None
            return None
            raise logger.error('Found too many data rows in out info file: %s' % file_obj.filename)

        file_obj.data = numpy.zeros((1, len(new_lines[0])), dtype=numpy.chararray)
       
        file_obj.data[0, :] = numpy.array(new_lines[0])
        file_obj.dims = file_obj.data.shape

    if file_obj.data == None:
        logger.error('Invalid outinfo file: %s for run dir: %s' % (file_obj.filename, curr_dir))
        return None

    try:
        return extract_file_value(ts_obj, file_obj, curr_dir, -1, column_desc, data_label)
    except LookupError as e:
        logger.error('Problem reading out info file, %s' % e)
        return None
        
def apply_psurf_to_pressures(pressures, psurf):
    logger = logging.getLogger(LOGGING_NAME)
    
    new_n_levels = len(pressures)
    for lev_idx in range(1, len(pressures)):
        if psurf < pressures[lev_idx]:

            psurf_diff = psurf - pressures[lev_idx-1]
            level_frac = (pressures[lev_idx]-pressures[lev_idx-1]) * NEW_LEVEL_FRAC
          
            if psurf_diff < level_frac:
                new_n_levels = lev_idx
            else:
                new_n_levels = lev_idx + 1

            break
                
    pressures = pressures[:new_n_levels]
    pressures[new_n_levels-1] = psurf

    return pressures

def calculate_aod_values(in_press, out_press, aer_profiles, name_prefix, use_sub_types=None):

    if type(aer_profiles) == numpy.ndarray and len(aerosol.shape) == 2:
        resampled_profiles = numpy.zeros((len(out_press), aer_profiles.shape[1]), dtype=float)
    else:
        resampled_profiles = numpy.zeros((len(out_press), len(aer_profiles)), dtype=float)
        
    for prof_idx in range(resampled_profiles.shape[1]):
        if type(aer_profiles) == numpy.ndarray and len(aerosol.shape) == 2:
            curr_profile = aer_profiles[:,prof_idx]
        else:
            curr_profile = numpy.array(aer_profiles[prof_idx])

        # Check for log values
        if numpy.any(curr_profile < 0):
            curr_profile = numpy.exp(curr_profile)
            
        resampled_profiles[:, prof_idx] = OCO_MathUtil.resample_profile(in_press, curr_profile, out_press)

    layer_aod, layer_press = OCO_MathUtil.total_aod( out_press, resampled_profiles, return_layers=True )

    if use_sub_types == None:
        use_sub_types = AOD_TYPE_ORDER

    aod_vals_ret = []
    aod_vals_chk = {}
    for type_name in AOD_TYPE_ORDER:
        type_range = AOD_SUB_TYPES[type_name]
        
        press_high = numpy.where(layer_press >= type_range[0])
        press_low  = numpy.where(layer_press <= type_range[1])
        press_range = list(set.intersection(set(press_high[0]), set(press_low[0])))

        calc_aod = numpy.sum(layer_aod[press_range])

        # Place variables with same name as type name for checking after loop
        aod_vals_chk[type_name] = calc_aod

        # Put calculated values in format needed for table
        if type_name in use_sub_types:
            aod_vals_ret.append( ('%s_%s' % (name_prefix, type_name), calc_aod) )

    check_eval = AOD_CHECK_EVAL.format(**aod_vals_chk)
    if not eval(check_eval):
        raise Exception('Check on aod values "%s" evaluated as: "%s" has failed.' % (AOD_CHECK_EVAL, check_eval))

    return aod_vals_ret
        
def extract_file_aod(ts_obj, file_obj, curr_dir, pressure_file, psurf_file, label_prefix, aerosol_re, use_sub_types=None):
    press_obj = ts_obj.get_cached_file(curr_dir, pressure_file, required=True)
    psurf_obj = ts_obj.get_cached_file(curr_dir, psurf_file)

    try:
        in_press_vals  = file_obj[PRESSURE_COLUMN][:,0]
    except LookupError:
        logger.error('For file "%s" LookupError: %s' % (file_obj.filename,traceback.format_exception_only(*sys.exc_info()[0:2])))
        return None

    if press_obj != None:
        out_press_vals = press_obj[PRESSURE_COLUMN][:,0]
    else:
        out_press_vals = in_press_vals

    try:
        if psurf_obj != None:
            out_press_vals = apply_psurf_to_pressures(out_press_vals, psurf_obj[PSURF_COLUMN][0,0])
    except LookupError:
        logger.error('For file "%s" LookupError: %s' % (press_obj.filename,traceback.format_exception_only(*sys.exc_info()[0:2])))
        return None

    if not hasattr(aerosol_re, '__iter__'):
        aerosol_re = ( aerosol_re, )

    aer_columns = []
    for col_name in file_obj.labels_lower:
        for curr_re in aerosol_re:
            if re.search(curr_re.lower(), col_name) and not re.search(PRESSURE_COLUMN.lower(), col_name):
                aer_columns.append( file_obj[col_name] )

    aod_data_vals = calculate_aod_values(in_press_vals, out_press_vals, aer_columns, label_prefix, use_sub_types)

    return aod_data_vals

def get_sv_names(ts_obj, curr_dir):
    logger = logging.getLogger(LOGGING_NAME)
    
    sv_obj = ts_obj.get_cached_file(curr_dir, SV_NAMES_BASE, as_strings=True)

    if sv_obj == None:
        logger.error('Could not load state vectors name file: %s for run dir: %s' % (SV_NAMES_BASE, curr_dir))
        return None

    ret_names = []
    for curr_name in sv_obj[SV_NAMES_COLUMN][:,0]:
        ret_names.append( str(curr_name).upper() )
        
    return ret_names

def get_sv_values(ts_obj, curr_dir, sv_obj, sv_item, sv_col, return_names=False):
    logger = logging.getLogger(LOGGING_NAME)
    
    sv_names = get_sv_names(ts_obj, curr_dir)

    if sv_names == None:
        logger.error('Could not retrieve statevector item: %s since statevector names file not available for run directory: %s' % (sv_item, curr_dir))
        return None

    col_values = sv_obj[sv_col][:,0]

    if len(sv_names) != len(col_values):
        logger.error('Statevector names size: %d and values size: %s do not match for run directory: %s' % (len(sv_names), len(col_values), curr_dir))
        return None
    
    sv_values = []
    used_names = []
    for curr_idx, curr_name in enumerate(sv_names):
        if re.search(sv_item, curr_name):
            sv_values.append( col_values[curr_idx] )
            used_names.append(curr_name)

    if return_names:
        return sv_values, used_names
    else:
        return sv_values
            
def extract_sv_item(ts_obj, file_obj, curr_dir, sv_item, item_idx, sv_col, data_label, report_missing=True, index_format=None):
    logger = logging.getLogger(LOGGING_NAME)
            
    item_values = get_sv_values(ts_obj, curr_dir, file_obj, sv_item, sv_col)

    if item_values == None or len(item_values) == 0:
        if report_missing:
            logger.error('Could not retrieve statevector item: %s for run directory: %s' % (sv_item, curr_dir))
        return None

    if item_idx == None:
        return_values = []
        if index_format != None:
            label_format = '%s'+index_format
        else:
            label_format = '%s_%d'
        for curr_idx, curr_item in enumerate(item_values):
            return_values.append( (label_format % (data_label, curr_idx+1), curr_item) )
        return return_values
    else:
        return [ (data_label, item_values[item_idx] ) ]

def extract_sv_aod(ts_obj, file_obj, curr_dir, pressure_file, sv_col, label_prefix, aerosol_re):
    logger = logging.getLogger(LOGGING_NAME)
    
    press_obj = ts_obj.get_cached_file(curr_dir, pressure_file, use_first=True)

    if press_obj == None:
        logger.error('Could not load pressure file: %s for run dir: %s' % (pressure_file, curr_dir))
        return None

    # If using out/atmosphere*.dat then retrieved psurf should already be in file, but just in case
    # or if you use a different file for pressure file here
    sv_psurf = get_sv_values(ts_obj, curr_dir, file_obj, PSURF_COLUMN, sv_col)

    try:
        if sv_psurf != None and len(sv_psurf) > 0:
            if sv_psurf[0] < VALID_PSURF_MIN:
                logger.error('Extracted invalid surface pressure value from statevector: %s. No statevector AOD computed for run dir: %s' % (sv_psurf, curr_dir))
                return None

            press_vals = apply_psurf_to_pressures(press_obj[PRESSURE_COLUMN][:,0], sv_psurf[0])
        else:
            press_vals = press_obj[PRESSURE_COLUMN][:,0]
    except LookupError:
        logger.error('For file "%s" LookupError: %s' % (press_obj.filename,traceback.format_exception_only(*sys.exc_info()[0:2])))
        return None

    if not hasattr(aerosol_re, '__iter__'):
        aerosol_re = ( aerosol_re, )

    aer_profiles = []
    for curr_re in aerosol_re:
        sv_values_ret = get_sv_values(ts_obj, curr_dir, file_obj, curr_re, sv_col, return_names=True)
        if sv_values_ret == None:
            continue
        
        aer_values, aer_names = sv_values_ret

        if len(aer_values) > len(press_vals):
            # First count how many types are listed, check how many of first index available
            num_types = 0
            for curr_name in aer_names:
                if re.search('_001$', curr_name):
                    num_types += 1

            if len(aer_values) != len(press_vals) * num_types:
                raise Exception('Calculated %d aerosol types for %d levels but only %d elements in available data for file: %s' % (num_types, len(press_vals), len(aer_values), file_obj.filename))

            aer_matrix = numpy.zeros((len(press_vals), num_types), dtype=float)

            sv_idx = 0
            for type_idx in range(num_types):
                for press_idx in range(len(press_vals)):
                    aer_matrix[press_idx, type_idx] = aer_values[sv_idx]
                    sv_idx += 1

            # Place aer_profiles list per type so values are parsed correctly
            for type_idx in range(num_types):
                aer_profiles.append( aer_matrix[:, type_idx] )
                
        elif len(aer_values) == len(press_vals):
            # Just one profile of the searched type
            aer_profiles.append( aer_values )
                   
        elif len(aer_values) > 0:
            logger.error('Unrecognized grouping of retrieval aerosols matching regular expression: "%s" for run directory: "%s". Number of aerosol values: %d, number of pressure levels: %d' % (curr_re, curr_dir, len(aer_values), len(press_vals)))

    return calculate_aod_values(press_vals, press_vals, aer_profiles, label_prefix)

def count_num_levels(ts_obj, file_obj, curr_dir, psurf_file, label_prefix):
    logger = logging.getLogger(LOGGING_NAME)
    
    psurf_obj = ts_obj.get_cached_file(curr_dir, psurf_file)

    try:
        in_press_vals  = file_obj[PRESSURE_COLUMN][:,0]
    except LookupError:
        logger.error('For file "%s" LookupError: %s' % (file_obj.filename,traceback.format_exception_only(*sys.exc_info()[0:2])))
        return None

    if psurf_obj != None:
        out_press_vals = apply_psurf_to_pressures(file_obj[PRESSURE_COLUMN][:,0], psurf_obj[PSURF_COLUMN][0,0])
    else:
        out_press_vals = file_obj[PRESSURE_COLUMN][:,0]

    return [ (label_prefix, len(out_press_vals)) ]

def compute_snr(ts_obj, file_obj, curr_dir, percent_range, num_highest, label_prefix):

    beg_idxs = file_obj.pixels
    end_idxs = file_obj.pixels[1:] + [file_obj.dims[0]]

    snr_data = []
    for band_idx, curr_beg, curr_end in zip(range(len(beg_idxs)), beg_idxs, end_idxs):       
        band_radiance = file_obj[RADIANCE_COLUMN][curr_beg:curr_end, 0]
        band_error    = file_obj[ERROR_COLUMN][curr_beg:curr_end, 0]

        # Get index for looking for SNR over from percentages supplied
        window_len = len(band_radiance)
        snr_range_beg = int(window_len*percent_range[0]/100.0)
        snr_range_end = int(window_len*percent_range[1]/100.0)

        band_snr = band_radiance[snr_range_beg:snr_range_end] / band_error[snr_range_beg:snr_range_end]
        
        # average N highest
        band_snr.sort()
        snr_mean_val = numpy.mean(band_snr[-num_highest:])

        snr_data.append( ('%s_%d' % (label_prefix, band_idx+1), snr_mean_val ) )

    return snr_data

def compute_interval_chi2(ts_obj, file_obj, curr_dir, meas_file, chi2_range, data_label):
    logger = logging.getLogger(LOGGING_NAME)
    
    meas_obj = ts_obj.get_cached_file(curr_dir, meas_file)

    if meas_obj == None:
        logger.error('Can not compute interval chi^2 since measured radiance file does not exist for run: %s' % curr_dir)
        return None

    try:
        wl_conv = file_obj[WAVELENGTH_COLUMN][:,0]
    except LookupError:
        logger.error('Could not retrieve column: %s from file: %s' % (WAVELENGTH_COLUMN, file_obj.filename))
        return None

    try:
        wl_meas = meas_obj[WAVELENGTH_COLUMN][:,0]
    except LookupError:
        logger.error('Could not retrieve column: %s from file: %s' % (WAVELENGTH_COLUMN, meas_obj.filename))
        return None

    conv_beg = bisect.bisect(wl_conv, chi2_range[0])
    conv_end = bisect.bisect(wl_conv, chi2_range[1])
    
    meas_beg = bisect.bisect(wl_meas, chi2_range[0])
    meas_end = bisect.bisect(wl_meas, chi2_range[1])

    residual = meas_obj[RADIANCE_COLUMN][meas_beg:meas_end, 0] - file_obj[RADIANCE_COLUMN][conv_beg:conv_end, 0]

    interval_chi2 = numpy.sum( numpy.power(residual / meas_obj[ERROR_COLUMN][meas_beg:meas_end, 0], 2) )

    return [ (data_label, interval_chi2) ]

SUMMARY_FILE_ID = 'Summary of L2 FP Runs'
# Dictionary is set up as follows:
# filename: ( ...extraction definition...)
# where extraction definition is a set of indexes and tuples
# if the item is an integer then it sets the insertion index for the columns defined afterwards until the next index.
# if the item is a tuple then it defines a function to call where the first item is the function name and the rest are arguments for the program
SUMMARY_FILE_CONFIG = \
                    { RESULT_FILE_BASE:           ( 10,
                                                    { 'as_strings': True },
                                                    (extract_result_value, -1, 3, 'YEAR'),
                                                    (extract_result_value, -1, 4, 'DAY'),
                                                    (extract_result_value, -1, 5, 'ZPD_TIME'), 
                                                    (extract_result_value, -1, 6, 'LATITUDE'), 
                                                    (extract_result_value, -1, 7, 'LONGITUDE'),

                                                    70,
                                                    (extract_result_value, -1, 8,  'SOLAR_ZENITH_1'),
                                                    (extract_result_value, -1, 18, 'SOLAR_ZENITH_2'),
                                                    (extract_result_value, -1, 28, 'SOLAR_ZENITH_3'),
                                                    (extract_result_value, -1, 10, 'OBS_ZENITH_1'),
                                                    (extract_result_value, -1, 20, 'OBS_ZENITH_2'),
                                                    (extract_result_value, -1, 30, 'OBS_ZENITH_3'),
                                                    ),

                      'true/x_target*.dat':        ( 20,
                                                     (extract_file_value, 0, 'TRUE_X_TARGET', 'REPORTED_XTARG_TRUE'),
                                                     ),
                      'true/xco2*.dat':            ( 20,
                                                     (extract_file_value, 0, 'TRUE_X_TARGET', 'REPORTED_XTARG_TRUE'),
                                                     ),
                      'out/control1/x_target.dat': ( 30,
                                                     (extract_file_value, 0, 'x_target', 'XTARG_RET'),
                                                     (extract_file_value, 0, 'a_priori', 'XTARG_RET_AP'),
                                                     (extract_file_value, 0, 'error',    'XTARG_RET_VAR'),
                                                     ),
                      'true/aerosol*.dat':         ( 40,
                                                     { 'use_first': True },
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_ALL',    '.*'     ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_ICE',    'ic_?.*' ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_WATER',  'wc_?.*' ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_AERO1',  'aero1'  ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_AERO2',  'aero2'  ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_0', 'clat_0', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_1', 'clat_1', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_2', 'clat_2', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_3', 'clat_3', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_4', 'clat_4', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_5', 'clat_5', ('TOTAL_AOD',) ),
                                                     (extract_file_aod, PRESSURE_TRUE_BASE, PSURF_TRUE_BASE, 'CALC_TRUE_CLAT_6', 'clat_6', ('TOTAL_AOD',) ),
                                                     ),
                      'true/aerosol_od*.dat':      ( 50,
                                                     (extract_file_value, 0, 'ICE_OD',     'REPORTED_ICE_AOD'),
                                                     (extract_file_value, 0, 'WATER_OD',   'REPORTED_WATER_AOD'),
                                                     (extract_file_value, 0, 'AEROSOL_OD', 'REPORTED_CONT_OCEANIC_AOD'),
                                                     ),

                      'out/aerosol_profiles.dat':  ( 60, # Use aerosol files written by L2 code instead of statevector values
                                                         # Since retrieval might be shape or profile type
                                                     (extract_file_aod, None, None, 'RET_ALL',    ('IC_?.*',
                                                                                                                'WC_?.*',
                                                                                                                'CONT', 'CLAT_0', 'AERO1',
                                                                                                                'OCEANIC','CLAT_1', 'AERO2')),
                                                     (extract_file_aod, None, None, 'RET_ICE',     'IC_?.*' ),
                                                     (extract_file_aod, None, None, 'RET_WATER',   'WC_?.*' ),
                                                     (extract_file_aod, None, None, 'RET_CONT',   ('CONT', 'CLAT_0', 'AERO1') ),
                                                     (extract_file_aod, None, None, 'RET_OCEANIC',('OCEANIC', 'CLAT_1', 'AERO2') ),

                                                     ),
                                                
                      STATEVECTOR_BASE:            ( #65, # Old method for 
                                                     #(extract_sv_aod, PRESSURE_RET_BASE, SV_RESULT_COL, 'RET_ALL',    ('IC_?.*',
                                                     #                                                                  'WC_?.*',
                                                     #                                                                  'CONT', 'CLAT_0', 'AERO1',
                                                     #                                                                  'OCEANIC','CLAT_1', 'AERO2')),
                                                     #(extract_sv_aod, PRESSURE_RET_BASE, SV_RESULT_COL, 'RET_ICE',     'IC_?.*' ),
                                                     #(extract_sv_aod, PRESSURE_RET_BASE, SV_RESULT_COL, 'RET_WATER',   'WC_?.*' ),
                                                     #(extract_sv_aod, PRESSURE_RET_BASE, SV_RESULT_COL, 'RET_CONT',   ('CONT', 'CLAT_0', 'AERO1') ),
                                                     #(extract_sv_aod, PRESSURE_RET_BASE, SV_RESULT_COL, 'RET_OCEANIC',('OCEANIC', 'CLAT_1', 'AERO2') ),
                                                     90,
                                                     (extract_sv_item, PSURF_COLUMN, 0, SV_APRIORI_COL, 'AP_PSURF'),
                                                     (extract_sv_item, PSURF_COLUMN, 0, SV_RESULT_COL,  'RET_PSURF'),
                                                     
                                                     110,
                                                     (extract_sv_item, 'WINDSPEED', 0, SV_RESULT_COL, 'RET_WINDSPEED', False),
                                                     
                                                     130,
                                                     (extract_sv_item, 'ALBEDO_001', 0, SV_RESULT_COL, 'RET_ALBEDO_1'),
                                                     (extract_sv_item, 'ALBEDO_002', 0, SV_RESULT_COL, 'RET_ALBEDO_2'),
                                                     (extract_sv_item, 'ALBEDO_003', 0, SV_RESULT_COL, 'RET_ALBEDO_3'),
                                                     200,
                                                     (extract_sv_item, 'EXP_01', None, SV_RESULT_COL, 'AEROSOL_EXP_01', False, '_%02d'),
                                                     (extract_sv_item, 'GAU_01', None, SV_RESULT_COL, 'AEROSOL_GAUSS_01', False, '_%02d'),
                                                     
                                                     ),
                      
                      PSURF_TRUE_BASE       :      ( 80,
                                                     (extract_file_value, 0, PSURF_COLUMN, 'TRUE_PSURF'),
                                                     ),
                      'true/windspeed*.dat':       ( 100,
                                                     (extract_file_value, 0, 'WINDSPEED', 'TRUE_WINDSPEED'),
                                                     ),
                      'true/albedo*.dat':          ( 120,
                                                     (extract_file_value, 0, 'ALBEDO_1', 'TRUE_ALBEDO_1'),
                                                     (extract_file_value, 0, 'ALBEDO_2', 'TRUE_ALBEDO_2'),
                                                     (extract_file_value, 0, 'ALBEDO_3', 'TRUE_ALBEDO_3'),
                                                     ),
                      ATMOSPHERE_TRUE_BASE:        ( 140,
                                                     (count_num_levels, PSURF_TRUE_BASE, 'NUM_TRUE_LEVELS'),
                                                     25,
                                                     (calculate_file_xco2, PSURF_TRUE_BASE, None, 'CALC_XTARG_TRUE'),
                                                     ),
                      PRESSURE_RET_BASE:           ( 160,
                                                     { 'use_first': True },
                                                     (count_num_levels, None, 'NUM_RET_LEVELS'),
                                                     # Removed since final atmosphere.dat does not have
                                                     # last state vector updates
                                                     #28,
                                                     #(calculate_file_xco2, None, None, 'CALC_XTARG_RET'),
                                                     ),
                      }

def noop(*args):
    pass

# Stats file
STATS_FILE_ID = 'Statistics of L2 FP Runs'
STATS_FILE_CONFIG = {  'out/aggregator/scalar.dat': ( 10,
                                                      (extract_header_value, 'abs_res_rms_o2',         'rms_abs_1'),
                                                      (extract_header_value, 'abs_res_rms_weak_co2',   'rms_abs_2'),
                                                      (extract_header_value, 'abs_res_rms_strong_co2', 'rms_abs_3'),
                                                      (extract_header_value, 'rel_res_rms_o2',         'rms_rel_1'),
                                                      (extract_header_value, 'rel_res_rms_weak_co2',   'rms_rel_2'),
                                                      (extract_header_value, 'rel_res_rms_strong_co2', 'rms_rel_3'),
                                                      30,
                                                      (extract_header_value, 'meas_chi2_norm_o2',         'chi2_1'),
                                                      (extract_header_value, 'meas_chi2_norm_weak_co2',   'chi2_2'),
                                                      (extract_header_value, 'meas_chi2_norm_strong_co2', 'chi2_3'),
                                                      (extract_header_value, 'tot_apriori_chi2',          'tot_apriori_chi2'), 
                                                      (extract_header_value, 'tot_measured_chi2',         'tot_measured_chi2'),
                                                      (extract_header_value, 'fc_tot_apriori_chi2',       'fc_tot_apriori_chi2'),
                                                      (extract_header_value, 'fc_tot_measured_chi2',      'fc_tot_measured_chi2'),
                                                      55,
                                                      (extract_header_value, 'iterations',                'iterations',),
                                                      (extract_header_value, 'divergences',               'divergences',),
                                                      ),
                       RAD_MEAS_BASE:               ( 20,
                                                      (compute_snr, (0,10),   2, 'beg_mean_snr'),
                                                      (compute_snr, (90,100), 2, 'end_mean_snr'),
                                                      ),
                       'out/rad_conv.dat':          ( 40,
                                                      (compute_interval_chi2, RAD_MEAS_BASE, ( 0.7625, 0.7655 ), 'sub_interval_chi2_o2'),
                                                      (compute_interval_chi2, RAD_MEAS_BASE, ( 1.602,  1.6075 ), 'sub_interval_chi2_weak_co2'),
                                                      (compute_interval_chi2, RAD_MEAS_BASE, ( 2.051,  2.056  ), 'sub_interval_chi2_strong_co2'),
                                                      ),
                       'out/out_info*.dat':         ( 50,
                                                      { 'as_strings': True },
                                                      (extract_outinfo_value, 'd_sigma_sq', 'd_sigma_sq'),
                                                      60,
                                                      (extract_outinfo_value, 'outcome',    'outcome'),
                                                      (extract_outinfo_value, 'conv_flag',  'conv_flag'),
                                                      (extract_outinfo_value, 'div_flag',   'div_flag'),
                                                      ),
                       }
                       
class testset_summary():

    def __init__(self, **addl_args):
        pass

    def load_existing_data(self, data_file):
        existing_obj = OCO_Matrix(data_file, as_strings=True)

        dict_data = existing_obj.as_dict(RUN_COLUMN_NAME)

        return existing_obj.labels, dict_data

    def convert_data_to_matrix(self, data_columns, data_dict):

        row_names = data_dict.keys()
        row_names.sort()

        matrix_data = numpy.zeros((len(row_names), len(data_columns)), dtype=numpy.chararray)
        matrix_data[:,:] = FILL_VALUE

        for row_idx, row_key in enumerate(row_names):
            for col_idx, col_name in enumerate(data_columns):

                data_row_dict = data_dict[row_key]
                if data_row_dict.has_key(col_name):
                    matrix_data[row_idx, col_idx] = str(data_row_dict[col_name])

        return matrix_data

    def reset_file_cache(self):
        self.file_cache = {}

    def get_cached_file(self, file_path, file_glob, required=False, use_first=False, **addl_args):
        if file_glob == None:
            return None
        
        if hasattr(file_glob, '__iter__'):
            found_files = []
            for curr_glob in file_glob:
                curr_files = glob.glob( os.path.join(file_path, curr_glob) )
                if len(curr_files) > 0:
                    found_files = curr_files
                    break
        else:
            found_files = glob.glob( os.path.join(file_path, file_glob) )

        if len(found_files) == 0:
            if required:
                raise OSError('Could not find at path: "%s" any files matching glob: "%s"' % (file_path, file_glob))
            else:
                return None
        elif len(found_files) > 1 and not use_first:
            raise OSError('Found too many files files at path: "%s" with glob: "%s", found: %s' % (file_path, file_glob, found_files))

        if self.file_cache.has_key(found_files[0]):
            file_obj = self.file_cache[ found_files[0] ]
        else:
            file_obj = OCO_Matrix( found_files[0], **addl_args )
            self.file_cache[ found_files[0] ] = file_obj

        return file_obj

    def extract_run_data(self, extract_config, run_dirs, run_names, output_filename, file_id=None, overwrite=False):
        logger = logging.getLogger(LOGGING_NAME)
        logger.info('')
        
        # Load existing data if possible and desired
        if not overwrite and os.path.exists(output_filename):
            logger.info('Loading existing data from: %s' % output_filename)
            column_names, data_dict = self.load_existing_data(output_filename)
        else:
            logger.info('Starting with fresh data dictionary')
            column_names = []
            data_dict = {}

        if RUN_COLUMN_NAME not in column_names:
            column_names.insert(0, RUN_COLUMN_NAME)

        # Record which names were encountered to determine if existing data
        # should remain in the file
        seen_names = []
            
        logger.info('Processing run directories:')
        for curr_idx, curr_name, curr_dir in zip(range(len(run_names)), run_names, run_dirs):
            if not os.path.exists(curr_dir):
                raise OSError('Run directory does not exist: %s' % curr_dir)

            # Mark this run as haven been seen and to be kept in the file
            seen_names.append(curr_name)

            if data_dict.has_key(curr_name):
                logger.info('%s skipped (%d of %d)' % (curr_name, curr_idx+1, len(run_names)))
                continue
            else:
                logger.info('%s <- %s (%d of %d)' % (curr_name, curr_dir, curr_idx+1, len(run_names)))
                data_dict[curr_name] = {}
                data_dict[curr_name][RUN_COLUMN_NAME] = curr_name
           
            self.reset_file_cache()
            insert_index = 0
            for extract_file, extract_items in extract_config.items():
                read_options = {}
                
                # Find read options
                for curr_item_desc in extract_items:
                    if type(curr_item_desc) == dict:
                        read_options = curr_item_desc

                extract_obj = self.get_cached_file(curr_dir, extract_file, **read_options)

                if extract_obj != None:
                    logger.debug('... %s' % extract_file)
                    for curr_item_desc in extract_items:
                        if type(curr_item_desc) == int:
                            insert_index = curr_item_desc
                        elif type(curr_item_desc) == dict:
                            # Ignore, already read in other loop
                            pass
                        elif type(curr_item_desc) == tuple:
                            try:
                                item_results = curr_item_desc[0](self, extract_obj, curr_dir, *curr_item_desc[1:])
                            except:
                                logger.error('Error calling routine: %s with values: %s' % (curr_item_desc[0], curr_item_desc[1:]))
                                logger.error(''.join(traceback.format_exception(*sys.exc_info(), limit=2)))
                                raise

                            # If nothing returned loop around
                            if item_results == None:
                                continue
                            
                            for item_name, item_value in item_results:
                                if item_name not in column_names:
                                    # Grow so insert indexes are meaningful 
                                    if insert_index > len(column_names)-1:
                                        column_names.extend([None] * (insert_index - len(column_names)+1))

                                    # If column name index is empty place item name there otherwise
                                    # after it
                                    if column_names[insert_index] == None:
                                        column_names[insert_index] = item_name
                                    else:
                                        column_names.insert(insert_index+1, item_name)
                                    
                                data_dict[curr_name][item_name] = item_value

                                # Increment so that items from current item results go in order
                                insert_index += 1

                        else:
                            raise Exception('Unknown item type: %s for extract file: %s' % (type(curr_item_desc), extract_file))
                else:
                    logger.debug('... %s (skipped, not present)' % extract_file)

            # Check validity of computed values
            for check_items, check_str in PER_RUN_CHECKS:
                
                check_count = 0
                for curr_item in check_items:
                    if curr_item in data_dict[curr_name].keys():
                        check_count += 1

                # Ignore if none of the required items are not in the dictionary for item
                if check_count == 0:
                    continue
                elif check_count != len(check_items):
                    err_msg = 'Not all required check items: "%s" present in dictionary for run_dir: %s' % (check_items, curr_dir)
                    logger.error(err_msg)
                    raise Exception(err_msg)
                else:
                    check_eval = check_str.format(**data_dict[curr_name])
                    if not eval( check_eval ):
                        err_msg = 'Check failed for run dir: %s. Check string: "%s" evaluated as: "%s"' % (curr_dir, check_str, check_eval)
                        logger.error(err_msg)
                        raise Exception(err_msg)

        # Cleanup empty values from column names:
        while True:
            try:
                column_names.remove(None)
            except ValueError:
                break

        # Remove any data that was not in the list of items
        # we were told to process, in case old data is present
        # in a file where the run dir list has been updated
        for data_name in data_dict.keys():
            if data_name not in seen_names:
                del data_dict[data_name]
                
        # Convert extracted data to matrix and write
        file_obj = OCO_Matrix()
        if file_id != None:
            file_obj.file_id = file_id
        file_obj.data = self.convert_data_to_matrix(column_names, data_dict)
        
        file_obj.labels = column_names
        file_obj.format_col = [True] * len(column_names)
        file_obj.format_col[ column_names.index(RUN_COLUMN_NAME) ] = False

        logger.info('Writing: %s' % output_filename)
        file_obj.write(output_filename, auto_size_cols=True)

    def testset_summary(self, run_dirs, run_names, summ_filename, overwrite=False):
        self.extract_run_data(SUMMARY_FILE_CONFIG, run_dirs, run_names, summ_filename, SUMMARY_FILE_ID, overwrite)

    def testset_stats(self, run_dirs, run_names, stat_filename, overwrite=False):
        self.extract_run_data(STATS_FILE_CONFIG, run_dirs, run_names, stat_filename, STATS_FILE_ID, overwrite)

    def make_summary_output(self, run_dirs, summ_filename, stat_filename, overwrite=False):
        logger = logging.getLogger(LOGGING_NAME)
        
        logger.debug('Extracting run names from run directories list')
        run_names = extract_run_names(run_dirs)

        self.testset_summary(run_dirs, run_names, summ_filename, overwrite=overwrite)
        self.testset_stats(run_dirs, run_names, stat_filename, overwrite=overwrite)

def standalone_main():

    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] [run_dir] [run_dir]...")

    parser.add_option( "-r", "--run_dir_file", dest="run_dir_file",
                       metavar="FILE",
                       help="file to read list of run directories from")

    parser.add_option( "--stat_file", dest="stat_filename",
                       metavar="FILE", default="./stats.dat",
                       help="name of stats table file other than default")

    parser.add_option( "--summ_file", dest="summ_filename",
                       metavar="FILE", default="./summary.dat",
                       help="name of summary table file other than default")

    parser.add_option( "--overwrite", dest="overwrite",
                       default=False,
                       action="store_true",
                       help="Overwrite existing summary and stats file instead of adding to existing contents"
                       )

    parser.add_option( "-v", "--verbose", dest="verbose",
                       default=False,
                       action="store_true",
                       help="Output more information to screen on processing"
                       )

    # Parse command line arguments
    (options, args) = parser.parse_args()

    # Initialize logging
    if options.verbose:
        L2_Log_Util.init_logging(logging.DEBUG)
    else:
        L2_Log_Util.init_logging(logging.INFO)

    # Gather list of run directories to gather information from
    run_dirs = []

    if (len(args) > 0):
        for arg_dir in args:
            run_dirs.append(arg_dir)

    if options.run_dir_file != None:
        if not os.path.exists(options.run_dir_file):
            parser.error("Run directory file '%s' does not exist" % options.run_dir_file)

        run_dir_fh = open(options.run_dir_file, 'r')
        for file_dir in run_dir_fh.readlines():
            run_dirs.append(file_dir.strip())
        run_dir_fh.close()

    # Sort items from run dir list
    run_dirs.sort()

    if len(run_dirs) == 0:
        parser.error('at least one run directory must be specified to be plotted')

    ts_summ = testset_summary()
    ts_summ.make_summary_output(run_dirs, options.summ_filename, options.stat_filename, overwrite=options.overwrite)

if __name__ == "__main__":
    standalone_main()
