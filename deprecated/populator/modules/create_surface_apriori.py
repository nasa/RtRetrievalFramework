import re
import math

import numpy

import ACOS_File
from OCO_Matrix import OCO_Matrix

from Generator_Utils import *

ALBEDO_COL_TMPL = 'ALBEDO_%d'
ALBEDO_FILE_TMPL = '%s/albedo_%s.dat'

BAND_CONTINUUM = { ACOS_File.OCO_INST_NAME:   ( range(18,32),
                              range(57,67),
                              [1,2] + range(37,39) ),
                   
                   ACOS_File.GOSAT_INST_NAME: ( range(2846, 2886) + range(3942, 3959),
                              range(4646, 4651) + range(4656, 4661),
                              ( range(2088, 2092), range(2273, 2277) + range(2484, 2488) ) )
                   }

BAND_DIODE_START = { ACOS_File.OCO_INST_NAME:   ( { 1016: 0 },
                                { 1016: 0 },
                                { 1016: 0 }
                                ),
                     ACOS_File.GOSAT_INST_NAME: ( { 6565:0, 1805:2356 },
                                { 8080:0, 3508:2256 },
                                { 6565:0, 2005:1755 }
                                ),
                     }

ALBEDO_CENTER_WAVELENGTHS = ( 0.77, 1.615, 2.06 )

MAX_ALBEDO_VALUE = 0.90

WINDSPEED_FG_BAND = 1
WINDSPEED_FG_PIXELS = range(57,68)
WINDSPEED_PARAM_NAMES = ('aod', 'beta', 'c0', 'c1', 'c2', 'c3')
WINDSPEED_FG_PARAMS = {
   (0.00, 0.05): dict(zip(WINDSPEED_PARAM_NAMES, (0.04, .004, -1.452, -4.577, -7.489, 6.357,))),
   (0.05, 0.15): dict(zip(WINDSPEED_PARAM_NAMES, (0.12, .016, -0.981, -4.925, -5.039, 4.560,))),
   (0.15, 0.25): dict(zip(WINDSPEED_PARAM_NAMES, (0.21, .016, -1.039, -3.434, -6.336, 4.958,))),
   }
WINDSPEED_FG_OFFSET = 0.6

WS_FG_A     = 9.054
WS_FG_B     = -14.42
WS_FG_ALPHA = 0.355
WS_FG_I0    = 1.8e21

# these need to be sequences
WINDSPEED_LABELS = 'WINDSPEED', 
WINDSPEED_UNITS  = 'm/s',

WINDSPEED_FILE_TMPL = '%s/windspeed_%s.dat'
WINDSPEED_F_FILE_TMPL = '%s/windspeed_f_%s.dat'

# For seeing if supplied aerosol value is a string or not
AEROSOL_FG_MATCH_RE = re.compile('0[0-9][0-9]')

# Information for extracting from an ASCII spectra file
ASCII_SPECTRA_RADIANCE_COLUMN = 'Radiance'
ASCII_SPECTRA_SZA_KEYWORD = 'sounding_solar_zenith'


def write_albedo_file(output_file, albedo_data, header_values=None):
   albedo_obj = OCO_Matrix()
   if header_values != None:
      albedo_obj.header.update(header_values)

   albedo_obj.header['center_wavelengths'] = ' '.join([str(wl) for wl in ALBEDO_CENTER_WAVELENGTHS])
   albedo_obj.labels = [ ALBEDO_COL_TMPL % (idx+1) for idx in range(albedo_data.shape[1]) ]
   albedo_obj.data = albedo_data
   albedo_obj.file_id = 'Surface albedo data'
   albedo_obj.write(output_file)

def create_albedo_apriori_from_radiance(radiance_data, sza_r, instrument_name, output_file):
   
   # Create apriori albedo based on continuum radiance
   # albedo=PI*I/(cos(SZA)*IO)
   # * 2.0 for oco for polarization sensitivity

   file_albedo_data = numpy.zeros((2, len(radiance_data)), dtype=float)
   orig_albedo_data = numpy.zeros(len(radiance_data), dtype=float)
   band_idx = 0
   for band_data in radiance_data:
      offset = BAND_DIODE_START[instrument_name][band_idx][len(band_data)]
      pixel_idxs = BAND_CONTINUUM[instrument_name][band_idx]
      
      max_ms_bnd = ACOS_File.MAX_MEAS_SIGNAL[instrument_name][band_idx]

      if type(pixel_idxs) is tuple:
         band_left, band_right = pixel_idxs
         band_left = [ idx-offset for idx in band_left ]
         band_right = [ idx-offset for idx in band_right ]
         cont_data  = [ max(band_data[band_left]), max(band_data[band_right]) ]
      else:
         pixel_idxs = [ idx-offset for idx in pixel_idxs ]
         cont_data  = band_data[pixel_idxs]
      albedo_val = math.pi * numpy.mean(cont_data) / (math.cos(sza_r) * max_ms_bnd)
      
      if instrument_name == ACOS_File.OCO_INST_NAME:
         albedo_val *= 2.0

      orig_albedo_data[band_idx] = albedo_val

      if albedo_val > 1.0:
         file_albedo_data[0, band_idx] = MAX_ALBEDO_VALUE
      else:
         file_albedo_data[0, band_idx] = albedo_val
      band_idx += 1

   header_values = {'calc_albedo' : ', '.join([ '%f' % aval for aval in orig_albedo_data])}
   write_albedo_file(output_file, file_albedo_data, header_values)

def create_sounding_albedo_apriori_from_radiance(sounding_id, radiance_data, sza_r, instrument_name, apriori_out_dir):
   albedo_filename = ALBEDO_FILE_TMPL % (apriori_out_dir, sounding_id)
   create_albedo_apriori_from_radiance(radiance_data, sza_r, instrument_name, albedo_filename)

def create_sounding_albedo_apriori_from_type_location(sounding_id, latitude, longitude, doy, scene_type, apriori_out_dir):
   file_albedo_data = numpy.zeros((2, len(ALBEDO_CENTER_WAVELENGTHS)), dtype=float)
   file_albedo_data[0,:] = Albedo_Table().get_albedo(doy, latitude, longitude, ALBEDO_CENTER_WAVELENGTHS, scene_type)

   write_albedo_file(sounding_id, apriori_out_dir, file_albedo_data)

def write_windspeed_file(sounding_id, apriori_out_dir, ws_value, header_values=None, windspeed_f=False):

   windspeed_obj = OCO_Matrix()

   windspeed_obj.data = numpy.zeros((1,1), dtype=float)
   windspeed_obj.data[:,:] = ws_value
   windspeed_obj.labels = WINDSPEED_LABELS
   windspeed_obj.units  = WINDSPEED_UNITS
   windspeed_obj.file_id = 'Windspeed for sounding: %s' % sounding_id

   if header_values != None:
      windspeed_obj.header.update(header_values)

   if windspeed_f:
      windspeed_obj.header['Windspeed_F'] = True
      windspeed_obj.write(WINDSPEED_F_FILE_TMPL % (apriori_out_dir, sounding_id))
   else:
      windspeed_obj.write(WINDSPEED_FILE_TMPL % (apriori_out_dir, sounding_id))
   
def create_windspeed_apriori_from_radiance(sounding_id, radiance_data, sza_r, aerosol_fg, apriori_out_dir, windspeed_f=False):

   # Convert aerosol first guess string if written in form where decimal is left out
   if isinstance(aerosol_fg, basestring) and re.match(AEROSOL_FG_MATCH_RE, aerosol_fg):
      aerosol_fg = float(aerosol_fg[0] + '.' + aerosol_fg[1:])

   aer_p = None # aersol parameters
   for aod_range, param_vals in WINDSPEED_FG_PARAMS.items():
      if aerosol_fg > aod_range[0] and aerosol_fg <= aod_range[1]:
         aer_p = param_vals
         break

   if aer_p == None:
      raise LookupError('Could not find valid windspeed first guess aerosol parameters for sounding id: %s with aerosol first guess: %s' % (sounding_id, aerosol_fg))
 
   i_wco2 = numpy.mean(radiance_data[WINDSPEED_FG_BAND][WINDSPEED_FG_PIXELS])
   mu = math.cos(sza_r)
   rscat = math.exp(aer_p['c0']+ aer_p['c1']*mu + aer_p['c2']*mu**2 + aer_p['c3'] * mu**3)
   rsurf = mu/(2.0*math.pi) * math.exp(WS_FG_A+WS_FG_B*mu**WS_FG_ALPHA) * math.exp(-2.0*aer_p['aod']/mu)
                        
   ws_f   = (i_wco2/WS_FG_I0 - rscat) / (aer_p['beta']*rscat + rsurf)
   ws_val = (1.0/ws_f - .003)/.00512 + WINDSPEED_FG_OFFSET

   ws_data = numpy.zeros((1,1), dtype=float)
   if windspeed_f:
      write_windspeed_file(sounding_id, apriori_out_dir, ws_f, windspeed_f=True) 
   else:
      write_windspeed_file(sounding_id, apriori_out_dir, ws_val, windspeed_f=False) 


def create_windspeed_apriori_from_ecmwf(sounding_id, ws_values, apriori_out_dir, windspeed_f=False):
   if windspeed_f:
      raise Exception('windspeed_f formulation not yet implemented')

   ws_u = ws_values[0][0]
   ws_v = ws_values[1][0]
   
   header_values = {}
   header_values['windspeed_u'] = ws_u
   header_values['windspeed_v'] = ws_v
   
   ws_val = numpy.sqrt( numpy.square(ws_u) + numpy.square(ws_v) )
   write_windspeed_file(sounding_id, apriori_out_dir, ws_val, header_values, windspeed_f=False) 

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):

    if len(moduleSections) > 1:
       raise Exception('Only one instance of %s allowed per FILE block' % os.path.basename(__file__))

    if source == destination:
       raise Exception('Will not write albedo into source spectra file. dest_filename must be defined')

    instrument_name = Apply_Template(moduleSections[0].Get_Keyword_Value('instrument_name'), valuesDict, mapDict=mapDict)

    if instrument_name == None or len(instrument_name) == 0:
        raise Exception('instrument_name keyword must be specified for module: %s' % os.path.basename(__file__))

    # Load spectra data source
    spec_file_obj = OCO_Matrix(source)

    radiances = []
    for pixel_range in spec_file_obj.pixel_ranges():
        radiances.append( spec_file_obj[ASCII_SPECTRA_RADIANCE_COLUMN][pixel_range[0]:pixel_range[1]] )

    # Get SZA value
    try:
       sza_r = math.radians( float(spec_file_obj.header[ASCII_SPECTRA_SZA_KEYWORD].split()[0]) )
    except KeyError:
       raise KeyError('Could not find header keyword %s in spectrum file: %s' % (ASCII_SPECTRA_SZA_KEYWORD, source))

    # Create apriori file from radiance data
    create_albedo_apriori_from_radiance(radiances, sza_r, instrument_name, destination)
