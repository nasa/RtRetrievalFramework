from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import os
import re
import sys
import math
import struct
import logging

my_list_type = None
try:
    # This is for python 2
    from types import ListType
    my_list_type = ListType
except ImportError:
    # This is for python 3
    my_list_type = list

from itertools import chain

import h5py
import numpy

from . import acos_file
from .math_util import *

SOUNDING_ID_COL_NAME = 'Sounding ID'

# Log albedo values
OS_ALBEDO_REF_POINTS = ( (0.755, 0.785),
                         (1.58, 1.65),
                         (2.03, 2.09), )
OS_ALBEDO_REF_POINTS = [ [ old_div(1e4, wl) for wl in band_wl ] for band_wl in OS_ALBEDO_REF_POINTS ]

L2_CENTER_WN = [ old_div(1e4, wvl) for wvl in ( 0.77, 1.615, 2.06 ) ]

EFF_ALBEDO_LOG_COLS = ( ('EffAlbedo_1A','EffAlbedo_1B',),
                        ('EffAlbedo_2A','EffAlbedo_2B',),
                        ('EffAlbedo_3A','EffAlbedo_3B',), )

# Log file filtering constants
AER_FILTER_COLS = ( 'Tau_water_%d',
                    'Tau_ice_%d',
                    'Tau_aerosol_%d', )
AER_FILTER_MAX = 0.3

logger = logging.getLogger()
      
class OrbitSimLogFile(object):
    def __init__(self, filename=None, l1b_file=None):

        self.filename = filename
        
        # Log files do not contain all the necessary
        # information for proper indexing, so need to
        # use the associated L1B file for help
        if l1b_file != None:
            self.l1b_obj = acos_file.L1B(l1b_file)
        else:
            self.l1b_obj = None
        
        self.comments = []
        self.column_names = []
        self.column_units = []
        self.data = None

        self.read()

        self.alb_indexes = []
        for spec_idx, band_col_names in enumerate(EFF_ALBEDO_LOG_COLS):
            self.alb_indexes.append([ self.get_column_index(nm) for nm in band_col_names])

        self.aer_indexes = []
        for spec_idx in range(len(self.alb_indexes)):
            self.aer_indexes.append([ self.get_column_index(aer_col % (spec_idx+1)) for aer_col in AER_FILTER_COLS ])

    def get_column_index(self, desired_col):

        if type(desired_col) is my_list_type:
            found_idx = []
            for curr_col_name in desired_col:
                found_idx.append( self.get_column_index(curr_col_name) )
            return found_idx
        else:
            column_name_srch = desired_col.strip()

            col_index = None
            if column_name_srch in self.column_names:
                col_index = self.column_names.index(column_name_srch)

            return col_index

    def get_sounding_indexes(self, sounding_id):
        if self.l1b_obj == None:
            raise LookupError("Can not return sounding indexes without a L1B file passed during object construction")

        return self.l1b_obj.get_sounding_indexes(sounding_id)

    def get_sounding_ids(self):
        if self.l1b_obj == None:
            raise LookupError("Can not return sounding ids without a L1B file passed during object construction")

        return numpy.ravel(self.l1b_obj.get_sounding_ids())

    def total_aod(self, frame_index, sounding_index, band_index):
        tot_aod = 0

        aer_indexes = self.aer_indexes[band_index]

        for curr_aer_idx in aer_indexes:
            tot_aod += self.data[frame_index, sounding_index, curr_aer_idx]
        return tot_aod

    def filter_retrievable(self, in_sounding_ids=None, return_indexes=False, max_aer=AER_FILTER_MAX):
        # Make sure is not a string
        max_aer = float(max_aer)

        if in_sounding_ids == None:
            in_sounding_ids = self.get_sounding_ids()

        out_soundings = []
        out_indexes = []
        for curr_sounding in in_sounding_ids:
            (frame_idx, snd_idx) = self.get_sounding_indexes(curr_sounding)
            FOOTPRINT_ID_COL_NAME = 'Frame'

            if self.total_aod(frame_idx, snd_idx, 0) < max_aer and \
                self.total_aod(frame_idx, snd_idx, 1) < max_aer and \
                self.total_aod(frame_idx, snd_idx, 2) < max_aer:
                out_soundings.append(curr_sounding)
                out_indexes.append(snd_idx)

        if return_indexes:
            return out_indexes
        else:
            return out_soundings

    def eff_albedo_coefficients(self, frame_index, sounding_index):
        eff_alb = numpy.zeros((len(self.alb_indexes), 2), dtype=float)
        for band_idx, band_alb_idxs in enumerate(self.alb_indexes):
            eff_alb[band_idx, :] = [ self.data[frame_index, sounding_index, ai] for ai in band_alb_idxs ]
        return eff_alb

    def l2_albedo_coefficients(self, frame_index, sounding_index, center_wn=None):
        eff_alb = self.eff_albedo_coefficients(frame_index, sounding_index)
        snd_alb_coeffs = numpy.zeros((eff_alb.shape[0], 2), dtype=float)

        for band_idx, (x_wn, c_wn) in enumerate(zip(OS_ALBEDO_REF_POINTS, L2_CENTER_WN)):
            coeff1 = old_div((eff_alb[band_idx, 1] - eff_alb[band_idx, 0]), (x_wn[1] - x_wn[0]))
            coeff0 = eff_alb[band_idx, 0] + (c_wn - x_wn[0]) * coeff1

            logger.debug("[%d] x0 = %e; x1 = %e" % (band_idx, x_wn[0], x_wn[1]))
            logger.debug("    y0 = %e; y1 = %e" % (eff_alb[band_idx, 0], eff_alb[band_idx, 1]))
            logger.debug("    c0 = %e; c1 = %e" % (coeff0, coeff1))

            snd_alb_coeffs[band_idx, :] = (coeff0, coeff1)

        return snd_alb_coeffs

    def read(self):
        logger = logging.getLogger(os.path.basename(__file__))
        
        col_name_detect = [ re.compile('#\s+Frame'), re.compile('#\s+Sounding') ]
        col_index_re = re.compile('(\s+\d)+')
        log_fo = open(self.filename, 'r')

        data_lines = []
        col_name_str = None
        col_unit_str = None
        col_idx_str = None
        for curr_line in log_fo.readlines():
            curr_line = curr_line.rstrip()

            #logger.debug('Parsing line: %s' % curr_line)

            found_col_name_line = False
            for curr_name_re in col_name_detect:
                if re.search(curr_name_re, curr_line):
                    found_col_name_line = True
                    break
            
            if found_col_name_line:
                logger.debug('Found column name line: %s' %curr_line)
                col_name_str = curr_line.lstrip('#')

            elif col_name_str != None and col_idx_str == None and re.search(col_index_re, curr_line):
                logger.debug('Found column index line: %s' % curr_line)
                col_idx_str = curr_line.lstrip('#')

            elif col_name_str != None and col_unit_str == None:
                logger.debug('Found column unit line: %s' % curr_line)
                col_unit_str = curr_line.lstrip('#')

            elif re.match('#', curr_line):
                self.comments = curr_line.lstrip('#')

            elif len(curr_line) > 0:
                data_lines.append(curr_line)
            
        log_fo.close()
        
        logger.debug("Read %d rows of data from file" % len(data_lines))

        if col_name_str == None or col_idx_str == None:
            raise IOError('Unable to find column label names in %s' % self.filename)

        column_lengths = []
        while len(col_idx_str) > 0:
            idx_match = re.search('\d+', col_idx_str)

            col_len = idx_match.end()+1
            if col_len < 0:
                raise IOError('Could not find next index in string "%s" for file %s' % (col_idx_str, self.filename))
            
            column_lengths.append( col_len )

            self.column_names.append( col_name_str[:col_len].strip() )
            self.column_units.append( col_unit_str[:col_len].strip() )

            col_idx_str = col_idx_str[col_len:]
            col_name_str = col_name_str[col_len:]
            col_unit_str = col_unit_str[col_len:]

        sound_id_col = self.get_column_index(SOUNDING_ID_COL_NAME)

        num_cols = len(self.column_names)
        num_rows = len(data_lines)

        self.data = numpy.zeros((num_rows, num_cols), dtype=float)

        soundings = []
        sounding_spec = {}
        for row_idx in range(num_rows):
            row_str = data_lines[row_idx]
            for col_idx in range(num_cols):
                col_len = column_lengths[col_idx]
                col_str = row_str[:col_len]
                row_str = row_str[col_len:]
                
                try:
                    self.data[row_idx, col_idx] = float(col_str)
                except ValueError:
                    raise ValueError('Could not convert: "%s" into a float at row: %d, column: %d in file: %s' % (col_str, row_idx, col_idx, self.filename))
          
                if col_idx == sound_id_col:
                    snd_id_val = int( self.data[row_idx, col_idx] )

                    if snd_id_val in sounding_spec:
                        sounding_spec[snd_id_val] += 1
                    else:
                        soundings.append(snd_id_val)
                        sounding_spec[snd_id_val] = 1

        if sound_id_col != None:
            num_rows = len(soundings)
            num_soundings = sounding_spec[soundings[0]]
        else:
            num_soundings = 1

        self.data.resize(num_rows, num_soundings, num_cols)

###

H2O_CONVERT_EPSILON = 0.621461
DRY_AIR_NAME = 'AIR_dry'

class SimulationSounding(object):

    def __init__(self, h5_obj, exposure_index):
        self.h5_obj = h5_obj
        self.exposure_index = exposure_index

    def close(self):
        self.h5_obj.close()

    def num_layer(self):
        return self.h5_obj["Simulation/Thermodynamic/num_layers"][self.exposure_index]

    def num_level(self):
        return self.num_layer() + 1

    def num_gas(self):
        return self.h5_obj["Simulation/Gas/num_species"][self.exposure_index]

    def gas_names(self):
        return [str(n) for n in self.h5_obj["Simulation/Gas/species_id"][self.exposure_index, :self.num_gas()]]

    def num_aerosol(self):
        return self.h5_obj["Simulation/Aerosol/num_species"][self.exposure_index]


    def aerosol_names(self):
        return [str(n) for n in self.h5_obj["Simulation/Aerosol/species_id"][self.exposure_index, :self.num_aerosol()]]

    def pressure(self):
        return self.h5_obj["Simulation/Thermodynamic/pressure_level"][self.exposure_index, 1:self.num_level()]
    
    def temperature(self):
        return self.h5_obj["Simulation/Thermodynamic/temperature_level"][self.exposure_index, 1:self.num_level()]


    def aerosol(self):
        aer_data = self.h5_obj["Simulation/Aerosol/species_density"][self.exposure_index, :self.num_aerosol(), :self.num_layer()]
        aer_dict = {}
        for aer_idx, aer_name in enumerate(self.aerosol_names()):
            aer_dict[aer_name] = aer_data[aer_idx, :]
        return aer_dict


    def gas(self):
        gas_data = self.h5_obj["Simulation/Gas/species_density"][self.exposure_index, :self.num_gas(), :self.num_layer()]
        gas_names = self.gas_names()
        dry_air_idx = gas_names.index(DRY_AIR_NAME)

        gas_dict = {}
        for gas_idx, gas_name in enumerate(gas_names):
            gas_val = gas_data[gas_idx, :]
            if gas_name.find("AIR") < 0:
                gas_val = old_div(gas_val, gas_data[dry_air_idx, :])

            gas_dict[gas_name] = gas_val
        return gas_dict

    def windspeed(self):
        return self.h5_obj["Simulation/Surface/wind_speed"][self.exposure_index]

class SimulationFile(h5py.File):

    NUM_EXPOSURES_DATASET = 'Simulation/Time/tai'
    
    def __init__(self, filename=None, l1b_reference=None):
        # Init super-class
        h5py.File.__init__(self, filename, 'r')
        
        self.l1b_ref = l1b_reference
        if self.l1b_ref != None:
            if not isinstance(self.l1b_ref, acos_file.L1B):
                raise TypeError("l1b_reference object must be of type acos_file.L1B")

            self.num_frames, self.num_soundings = self.l1b_ref.get_sounding_ids().shape[:2]

    def num_exposures(self):
        return self[self.NUM_EXPOSURES_DATASET].shape[0]

    def get_sounding_info(self, sounding_id):
        if self.l1b_ref == None:
            raise Exception("l1b_ref object was not passed in object construction")

        frame_idx, snd_idx = self.l1b_ref.get_sounding_indexes(sounding_id)
        exposure_index = frame_idx * self.num_soundings + snd_idx

        return self[exposure_index]

    def __getitem__(self, key):
        if isinstance(key, int):
            return Simulation_Sounding(self, key)
        else:
            return h5py.File.__getitem__(self, key)
