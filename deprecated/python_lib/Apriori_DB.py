import os
import time
import datetime

import numpy

from glob import glob

from OCO_Matrix import OCO_Matrix
import OCO_MathUtil

if not 'L2_DATASETS_PATH' in os.environ: raise Exception("L2_DATASETS_PATH not set, please see user's guide for setup instructions")

DEFAULT_TREND_DATABASE_FILE = os.path.join(os.environ['L2_DATASETS_PATH'], 'trends/co2_mm_gl.txt')
DEFAULT_LMDZ_APRIORI_PATH   = os.path.join(os.environ['L2_DATASETS_PATH'], 'apriori/co2')
DEFAULT_LMDZ_MEAN_FILE      = os.path.join(DEFAULT_LMDZ_APRIORI_PATH, 'lmdz_mean_values.dat')

# 2Describes how latitude files are binned
# Currently from -90 to 90 with 10 degree bins
LATITUDE_RANGE_BINS = numpy.array( range(-90, 91, 10) )

# Describes how apriori filenames are constructed
APRIORI_FILENAME_GLOB_TMPL = "*{num_levels}_{lat_bin_beg}_{lat_bin_end}_{month}_{surface_type}*.dat"

# Name of column to adjust in CO2 apriori files
CO2_APRIORI_CO2_COL = 'CO2'

# Default column to use from trend database file
DEFAULT_TREND_VALUES_COL = "trend"

# Scale to ppm
TREND_VALUE_SCALING = 1e-6

class Trend_DB_Data:
    def __init__(self, filename, values_column=DEFAULT_TREND_VALUES_COL):
        self.dates  = []
        self.values = []

        self.values_column = values_column

        self.read(filename)

    def add_value(self, value_dict):
        self.dates.append( float(value_dict['decimal']) )
        self.values.append( float(value_dict[self.values_column]) )

    def read(self, trend_db_file):
        """Load ESRL trend database ascii file"""

        self.filename = trend_db_file
        with open(trend_db_file) as trend_fobj:
            
            column_names = None
            for file_line in trend_fobj.readlines():
                line_parts = file_line.strip().split()

                if line_parts[0].find('#') == 0:
                    # Rest of line_parts could be column names
                    # Last comment before data should be column names
                    column_names = line_parts[1:]
                elif len(line_parts) > 0:
                    # If not a comment and not empty treat as data
                    line_dict = {}
                    for col_name, col_value in zip(column_names, line_parts):
                       line_dict[col_name] = col_value

                    self.add_value(line_dict)

    def get_value(self, time_struct):
        year, month, doy = (time_struct.tm_year, time_struct.tm_mon, time_struct.tm_yday)
        year_frac = float(doy) / float(time.strftime('%j', datetime.datetime(year, 12, 31).utctimetuple()))
        curr_dec_date = year + year_frac

        trend_value = OCO_MathUtil.linear_interpol(self.values, self.dates, [curr_dec_date], extrapolate=True, extrap_err=False)

        # Resulting trend value should be in first element of matrix
        return trend_value[0]

class Apriori_Mean_Data:
    """Class to load apriori mean information from file"""
    
    def __init__(self, filename):
        self.read(filename)

    def read(self, filename):
        self.filename = filename

        file_obj = OCO_Matrix(filename)

        self.global_mean = float(file_obj.header['global_ocean_mean'])

        self.months         = file_obj['Month']
        self.surface_values = file_obj['Surface_']
        
class CO2:
    def __init__(self, apriori_db_path=DEFAULT_LMDZ_APRIORI_PATH, apriori_mean_file=DEFAULT_LMDZ_MEAN_FILE, trend_db_file=DEFAULT_TREND_DATABASE_FILE, alternative_global_mean=None):
        self.apriori_db_path = apriori_db_path
        
        self.trend_db = Trend_DB_Data(trend_db_file)
        self.apriori_db_mean = Apriori_Mean_Data(apriori_mean_file)
        self.alternative_global_mean = alternative_global_mean

    def get_latitude_bin_bounds(self, latitude):
       """Get beginning and end latitudes for use in apriori filenames"""

       where_bin_beg = numpy.where(LATITUDE_RANGE_BINS >= latitude)
       try:
           lat_bin_beg = LATITUDE_RANGE_BINS[where_bin_beg[0][0]-1]
           lat_bin_end = LATITUDE_RANGE_BINS[where_bin_beg[0][0]]
       except IndexError:
           raise ValueError('Picked out of range bin for latitude: %f' % (latitude))

       return lat_bin_beg, lat_bin_end

    def get_apriori_offset(self, time_struct, debug_values={}):
       """
Get offset to be added to apriori co2 file according to this model

a_priori (year, month) =
   LMDZ profile(latitude band, month, levels)
   + ESRL surface co2 global mean trend value (year, month)
   - LMDZ surface co2 global ocean mean

       """

       if self.alternative_global_mean != None:
           current_global_mean = self.alternative_global_mean
       else:
           current_global_mean = self.trend_db.get_value(time_struct)

       offset = (current_global_mean - self.apriori_db_mean.global_mean)*TREND_VALUE_SCALING

       debug_values['current_global_mean']  = current_global_mean
       debug_values['database_global_mean'] = self.apriori_db_mean.global_mean
       debug_values['co2_offset']           = offset

       return offset
       
    def adjust_file_for_trend(self, src_apriori_file, time_struct, dst_apriori_file=None):

       if dst_apriori_file == None:
           dst_apriori_file = src_apriori_file

       if type(dst_apriori_file) is str and os.path.realpath(os.path.dirname(dst_apriori_file)) == os.path.realpath(self.apriori_db_path):
           raise IOError('Can not modify apriori file as located in database path, it must be copied first')

       apriori_obj = OCO_Matrix(src_apriori_file)
       co2_col_idx = apriori_obj.labels.index(apriori_obj.find_labels(CO2_APRIORI_CO2_COL)[0])

       co2_offset = self.get_apriori_offset(time_struct, debug_values=apriori_obj.header)
       apriori_obj.data[:, co2_col_idx] += co2_offset
      
       apriori_obj.header['co2_offset']  = co2_offset

       apriori_obj.write(dst_apriori_file)

    def find_apriori_files(self, latitude, time_struct, surf_type=None, num_levels=None):

       lat_bin_beg, lat_bin_end = self.get_latitude_bin_bounds(latitude)

       # Set the real default values
       if surf_type == None:
           surf_type = '*'

       if num_levels == None:
           num_levels= '*'

       search_format_dict = { 'num_levels':   num_levels,
                              'lat_bin_beg':  lat_bin_beg,
                              'lat_bin_end':  lat_bin_end,
                              'month':        time_struct.tm_mon,
                              'surface_type': surf_type,
                              }
       apriori_search_glob = APRIORI_FILENAME_GLOB_TMPL.format(**search_format_dict)

       # Allow apriori database path to contain multiple directories
       if hasattr(self.apriori_db_path, '__iter__'):
           data_source_dir = self.apriori_db_path
       else:
           data_source_dir = ( self.apriori_db_path, )

       matching_files = []
       for curr_source_dir in data_source_dir:
          if not os.path.exists(curr_source_dir) or not os.path.isdir(curr_source_dir):
              err_msg = 'data_source_dir: %s does not exist or is not a directory' % curr_source_dir
              raise IOError(err_msg)

          file_srch_glob = os.path.join(curr_source_dir, apriori_search_glob)
          found_files = glob(file_srch_glob)

          matching_files += found_files


       return found_files
