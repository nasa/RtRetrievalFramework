import os
import re
import glob
import shutil
import time
import datetime
import bisect
import math
import numpy
from contextlib import contextmanager, nested
from types import DictType
import itertools
import logging
import traceback

import h5py

from Generator_Utils import *
from OCO_TextUtils import evaluate_bool_str

import ACOS_File
import Apriori_DB

from load_co2_apriori import find_and_adjust_co2_apriori
from load_disp_apriori import create_scene_dispersion_file
from create_surface_apriori import create_sounding_albedo_apriori_from_radiance, create_windspeed_apriori_from_radiance, create_windspeed_apriori_from_ecmwf

# How to name CO2 apriori file
CO2_APRIORI_FILE_TMPL = 'co2_apriori_%s.dat'

# How to name dispersion output file
DISPERSION_FILE_TMPL = 'dispersion_%s.dat'

LAMBERTIAN = 'lambertian'
COXMUNK = 'coxmunk'

# Mapping of surface type to brdf method
APRIORI_SURF_MASK_PICKS = {
   'land'  : LAMBERTIAN,
   'ocean' : COXMUNK,
   }

WINDSPEED_U_DATASET = '/ecmwf/windspeed_u'
WINDSPEED_V_DATASET = '/ecmwf/windspeed_v'

@contextmanager
def write_if_def(filename, verbose_name=None):
   logger = logging.getLogger(os.path.basename(__file__))

   # __enter__
   if verbose_name != None:
      logger.debug('Opening %s file for writing: %s' % (verbose_name, filename))
      
   if filename != None and len(filename) != 0:
      file_obj = open(filename, 'w')
   else:
      file_obj = None

   # do some work
   yield file_obj

   # __exit__
   if file_obj != None:
      file_obj.close()

def create_static_apriori_files(sounding_id, latitude, longitude, time_struct, num_levels, surf_type, apriori_out_dir, co2_apriori_db):

   # Return the name of the file we found
   co2_apriori_file = os.path.join(apriori_out_dir, CO2_APRIORI_FILE_TMPL % sounding_id)
   source_file, dest_file = find_and_adjust_co2_apriori(co2_apriori_db, co2_apriori_file, latitude, time_struct, surf_type, num_levels)
   
   return os.path.basename(dest_file)

def load_sounding_apriori(sounding_id, l1b_obj, ecmwf_obj, num_levels, apriori_out_dir, co2_apriori_db, apriori_base_id_obj=None, brdf_type_id_obj=None, gain_type_id_obj=None, windspeed_f=False, l1b_build_id=None, use_brdf_type=None):

   if l1b_obj.instrument_name == ACOS_File.GOSAT_INST_NAME and type(sounding_id) is str and sounding_id[-1] not in ACOS_File.GOSAT_POL_ORDER:
      average_name = 'Polarization'
   else:
      average_name = None

   # Get radiances per band
   radiances   = l1b_obj.get_radiance_data(sounding_id, average=average_name)

   # These return per band, so just get first one
   sounding_indexes = l1b_obj.get_sounding_indexes(sounding_id)
   
   latitude    = l1b_obj.get_sounding_info('latitude', sounding_id, average=average_name)[0]
   longitude   = l1b_obj.get_sounding_info('longitude', sounding_id, average=average_name)[0]
   time_struct = l1b_obj.get_sounding_time(sounding_id, average=average_name)[0]
   sza_r       = math.radians(l1b_obj.get_sounding_info('solar_zenith', sounding_id, average=average_name)[0])
   saz_r       = math.radians(l1b_obj.get_sounding_info('solar_azimuth', sounding_id, average=average_name)[0])
   apriori_surf_type = l1b_obj.get_surface_grouping(sounding_id)

   if l1b_obj.instrument_name == ACOS_File.GOSAT_INST_NAME and gain_type_id_obj != None:
      gain_code = l1b_obj.get_sounding_info('gain', sounding_id)
      print >>gain_type_id_obj, "'%s' '%s'" % (sounding_id, gain_code)

   if use_brdf_type != None and len(use_brdf_type) > 0:
      brdf_type = use_brdf_type
   else:
      # Get the mapping of the apriori surface type to the albedo to be used
      brdf_type = APRIORI_SURF_MASK_PICKS[ apriori_surf_type ]

   if brdf_type == None:
      raise LookupError('Could not determine brdf type from surface mask')

   if brdf_type_id_obj != None:
      print >>brdf_type_id_obj, "%s %s" % (sounding_id, brdf_type)

   apriori_map_value = create_static_apriori_files(sounding_id, latitude, longitude, time_struct, num_levels, apriori_surf_type, apriori_out_dir, co2_apriori_db)

   if apriori_base_id_obj != None:
      print >>apriori_base_id_obj, "%s %s" % (sounding_id, apriori_map_value)

   if brdf_type == LAMBERTIAN:
      # Create albedo apriori from radiances
      create_sounding_albedo_apriori_from_radiance(sounding_id, radiances, sza_r, l1b_obj.instrument_name, apriori_out_dir)
   elif brdf_type == COXMUNK:
      # Create windspeed file
     ws_u_val = ecmwf_obj[WINDSPEED_U_DATASET][sounding_indexes[0],:,sounding_indexes[1]]
     ws_v_val = ecmwf_obj[WINDSPEED_V_DATASET][sounding_indexes[0],:,sounding_indexes[1]]
     if len(ws_u_val.shape) > 1:
        ws_u_val = numpy.average(ws_u_val, 1)
        ws_v_val = numpy.average(ws_v_val, 1)
     ecmwf_ws_values = (ws_u_val, ws_v_val)

     create_windspeed_apriori_from_ecmwf(sounding_id, ecmwf_ws_values, apriori_out_dir, windspeed_f)
         
   else:
      raise Exception("Unknown brdf_type: %s for sounding %s" % (brdf_type, sounding_id))

   if l1b_obj.instrument_name == ACOS_File.GOSAT_INST_NAME:
      dispersion_coefs = l1b_obj.get_sounding_info('dispersion', sounding_id, average=average_name)

      dispersion_ap_file = os.path.join(apriori_out_dir, DISPERSION_FILE_TMPL % sounding_id)

      create_scene_dispersion_file(sounding_id, latitude, sza_r, saz_r, time_struct, radiances[0], dispersion_coefs, dispersion_ap_file, l1b_build_id)

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
   logger = logging.getLogger(os.path.basename(__file__))
   
   if len(moduleSections) > 1:
      raise RuntimeError( 'Only one extraction block allowed')

   apriori_dir      = Apply_Template(moduleSections[0].Get_Keyword_Value('apriori_dir'), valuesDict, mapDict=mapDict)
   num_levels       = int(Apply_Template(moduleSections[0].Get_Keyword_Value('num_levels'), valuesDict, mapDict=mapDict))
   l1b_file         = Apply_Template(moduleSections[0].Get_Keyword_Value('l1b_file'), valuesDict, mapDict=mapDict)
   ecmwf_file       = Apply_Template(moduleSections[0].Get_Keyword_Value('ecmwf_file'), valuesDict, mapDict=mapDict)
   sounding_id_file = Apply_Template(moduleSections[0].Get_Keyword_Value('sounding_id_file'), valuesDict, mapDict=mapDict)
   sounding_id_sect = Apply_Template(moduleSections[0].Get_Keyword_Value('sounding_id_sect'), valuesDict, mapDict=mapDict)

   # Force all soundings to use a specific brdf type
   use_brdf_type    = Apply_Template(moduleSections[0].Get_Keyword_Value('use_brdf_type'), valuesDict, mapDict=mapDict)

   apriori_base_id_map = Apply_Template(moduleSections[0].Get_Keyword_Value('apriori_base_id_map'), valuesDict, mapDict=mapDict)
   brdf_type_id_map    = Apply_Template(moduleSections[0].Get_Keyword_Value('brdf_type_id_map'), valuesDict, mapDict=mapDict)
   gain_type_id_map    = Apply_Template(moduleSections[0].Get_Keyword_Value('gain_type_id_map'), valuesDict, mapDict=mapDict)
   windspeed_f         = evaluate_bool_str(Apply_Template(moduleSections[0].Get_Keyword_Value('windspeed_f'), valuesDict, mapDict=mapDict), default=False)

   global_mean_xco2    = Apply_Template(moduleSections[0].Get_Keyword_Value('global_mean_xco2'), valuesDict, mapDict=mapDict)

   # Check arguments
   if not os.path.isdir(apriori_dir):
      raise IOError('destFilename: %s must be a directory' % apriori_dir)

   for required_keyword in ('num_levels', 'l1b_file', 'sounding_id_file', 'sounding_id_sect'):
      if eval(required_keyword) == None:
       raise ValueError('%s keyword must be defined' % required_keyword)

   for file_check in ('l1b_file', 'sounding_id_file'):
      if not os.path.exists(eval(file_check)):
         raise IOError('%s: %s does not exist' % (file_check, eval(file_check)))

   # Load data from external files needed to only be read once
   sounding_id_list = Read_Id_List_File(sounding_id_file, sounding_id_sect, valuesDict=valuesDict, mapDict=mapDict)

   if global_mean_xco2 != None:
      if len(global_mean_xco2) == 0:
         global_mean_xco2 = None
      else:
         global_mean_xco2 = float(global_mean_xco2)

   co2_apriori_db = Apriori_DB.CO2(alternative_global_mean=global_mean_xco2)

   # Operate in context of input files
   with nested(ACOS_File.L1B(l1b_file), h5py.File(ecmwf_file,'r')) as (l1b_obj, ecmwf_obj):
       logger.debug('Loading apriori data based on %d soundings from l1b file: %s' % (len(sounding_id_list), l1b_file))

       # Retieve these here since they are done with H5Dump and should not be
       # done for each sounding
       l1b_build_id = l1b_obj.get_build_id()

       # Write per sounding files and maps in context of output map files
       with nested(write_if_def(apriori_base_id_map, 'apriori base name map'), write_if_def(brdf_type_id_map, 'brdf type map'), write_if_def(gain_type_id_map, 'gain code map'),) as (apriori_base_id_obj, brdf_type_id_obj, gain_type_id_obj):


           sys.stdout.write('%s: Creating per sounding data: ' % os.path.basename(__file__))
           for sounding_id in sounding_id_list:

              load_sounding_apriori(sounding_id, l1b_obj, ecmwf_obj, num_levels, apriori_dir, co2_apriori_db, apriori_base_id_obj, brdf_type_id_obj, gain_type_id_obj, windspeed_f, l1b_build_id=l1b_build_id, use_brdf_type=use_brdf_type)

              # Progress marks
              sys.stdout.write('.')
              sys.stdout.flush()

           # Newline for progress marks
           sys.stdout.write('\n')
