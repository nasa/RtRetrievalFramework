import os

from scipy import constants

import ACOS_File
import OCO_TextUtils
from OCO_Matrix import OCO_Matrix

from Generator_Utils import *

# Used for dispersion calculation
ABAND_SOLAR_LINE = 12985.16325e0 # line position of strong, isolated solar line
ABAND_ILS_OFFSET = { 'P': 0.203e0, 'S': 0.203e0 } # this takes into account the offset of the GOSAT band 1 ILS.

# Averaged offset when not using P or S seperately
ABAND_ILS_OFFSET['averaged'] = numpy.average(ABAND_ILS_OFFSET['P'], ABAND_ILS_OFFSET['S'])

ABAND_DISP_SEARCH_RANGE = [12983.0,12988.0] # search range

DISPERSION_OFFSET_SCALING = [ 1, 0.48, 0.38 ]
DISPERSION_ASCII_COLUMN_IDENT = 'DISP_'
DISPERSION_LABELS = [ 'PARAMETER' ] + [ '%s%d' % (DISPERSION_ASCII_COLUMN_IDENT, disp_idx) for disp_idx in range(1,4) ]

# For solar doppler correction
CON = 6.73951496e-03 # parameter for calculating geocentric from geodetic latitude
EARTH_ROT_SPD = 2.0*math.pi/86164.09054e0 # earth angular rotation frequency [1/sec]
RADIUS = 6378.137e0 # Earth equatorial radius, in km

def create_scene_dispersion_file(sounding_id, latitude, sza_r, saz_r, time_struct, aband_data, dispersion_coefs, apriori_out_file, l1b_build_id=None, index_scheme=None):
   """
   Purpose: Given a single frame of gosat data, estimate the spectral shift in band 1P.
            This is accomplished by using the strong solar line at 12985.163 wavenumbers.
            This is in the a-band continuum and is really the only local feature there.
   
            Note: This routine fits for a global shift, which INCLUDES the instrument
            doppler shift.  If this is not desired, the user must subtract off
            the instrument doppler shift.
   
   Original Author: Chris Odell
   Converted from IDL
   """

   if type(sounding_id) is str:
      pol_name = sounding_id[-1]
      if pol_name not in ACOS_File.GOSAT_POL_ORDER:
         # If not P or S then use averaged offset
         pol_name = 'averaged'
   else:
      raise Exception('Can not determine polarization name from non string ;Asounding id: %s' % sounding_id)

   wn_s0 = ABAND_SOLAR_LINE
   wn_s0 += ABAND_ILS_OFFSET[pol_name] # account for different ILS offset P vs S

   # Calculate the Speed of Earth center towards sun
   frac = (time_struct.tm_hour + time_struct.tm_min/60. + time_struct.tm_sec/3600.)/24. # fraction of a day
   doy = time_struct.tm_yday + frac + 0.5 # why 0.5? well cause L2 code does that
   V_cen = 497.2 * math.sin((doy-4.1)/365.25*2*math.pi) # simple approximate model, good to ~ 10 m/s.

   # Calculate the rotational component of the solar doppler velocity
   geo_lat = math.atan(math.tan(math.radians(latitude))/(1.+CON))
   Rloc = 1000.0 * RADIUS/math.sqrt(1.+ CON*math.sin(geo_lat)**2)
   V_rot = -EARTH_ROT_SPD * Rloc * math.sin(sza_r) * math.cos(saz_r-math.radians(90.)) * math.cos(geo_lat)

   solar_doppler_velocity = V_cen + V_rot
  
   # Modify position of strong solar line to include solar doppler shift
   wn_s = wn_s0 * (1.0e0 - solar_doppler_velocity/constants.c) # takes into account the solar doppler shift

   # Not used but could be good for debugging?
   #solar_shift = wn_s - wn_s0

   # Create aband dispersion array, 1 indexed
   abo2_dcoef_0, abo2_dcoef_1 = dispersion_coefs[0, :]
   aband_disp = numpy.arange(1, len(aband_data)+1)*abo2_dcoef_1 + abo2_dcoef_0

   # (2) CONSTRUCT THE X, Y FUNCTION TO FIT
   #
   # Use intersection of the points above and below range threshold
   # Make sure to sort the resulting indexes since set.intersection
   # does not guarantee anything about ordering
   w_1 = numpy.where(aband_disp >= ABAND_DISP_SEARCH_RANGE[0])
   w_2 = numpy.where(aband_disp <= ABAND_DISP_SEARCH_RANGE[1])
   w = sorted(list(set.intersection(set(w_1[0]), set(w_2[0]))))

   if len(w) == 0:
      raise ValueError('Could not find any points in the range %s in the aband dispersion with range %s for sounding %s' % (ABAND_DISP_SEARCH_RANGE, (min(aband_disp), max(aband_disp)), sounding_id))

   x = aband_disp[w] - wn_s

   mxmeas = max(aband_data[w]) # maximum value
   m2 = max(mxmeas - aband_data[w])
   y = (mxmeas-aband_data[w])/m2 # this should look like a gaussian

   w = numpy.where(y < 0.1)
   y = y - numpy.mean(y[w]) # subtract off the "continuum"
   pos = numpy.argmax(y) # pos = index of the maximum value of y
   m = len(y)
   fg = -x[pos]  # first-guess value of spectral shift

   # (3) PERFORM A SIMPLE FIT ASSUMING A PERFECTLY LINEAR MODEL
   #     MODEL IS A GAUSSIAN WITH SIGMA=0.2 CM^-1
   #     FIT PARAMETERS ARE AMPLITUDE AND CENTER OF THE GAUSSIAN
   K = numpy.zeros((m, 2), dtype=float) # will hold jacobian
   sig2 = 0.2e0**2
   ygauss = numpy.exp(-(x+fg)**2/ (2.0*sig2))
   K[:, 0] = ygauss
   K[:, 1] = -ygauss / sig2 * (x+fg)
   Kt = numpy.transpose(K)

   # Note that KtK is a 2x2 symmetric matrix and can be inverted analytically if desired.
   # form matrix (Kt K)^{-1} Kt
   delta = numpy.dot(numpy.linalg.inv(numpy.dot(Kt, K)), (numpy.dot(Kt, (y-ygauss))))
   fit = [1.0, fg] + delta

   wn_shift = fit[1]

   # Debugging to make sure fitting works
   if (0):
      from matplotlib import pyplot
      print sounding_id, aband_disp[0], wn_shift, aband_disp[0]-aband_disp[1]

      pyplot.cla()
      pyplot.plot(x, y) # measured data
      yfit = fit[0]* numpy.exp(-(x+fit[1])**2/ (2.0*sig2))
      pyplot.plot(x,yfit) # fit modeled data

      pyplot.legend(('Measured Solar Line', 'Fitted Solar Line'), loc='lower left')
      
      pyplot.savefig(os.path.join(os.path.basename(apriori_out_file), 'disp_fit_%s.png' % sounding_id), format='png')


   disp_obj = OCO_Matrix()

   disp_obj.data = numpy.zeros((dispersion_coefs.shape[1], dispersion_coefs.shape[0]+1), dtype=float)

   # First column is index of parameter
   disp_obj.data[:, 0] = numpy.arange(1,dispersion_coefs.shape[1]+1)

   # Create dispersion file scaling non aband values according to this method:
   # Where below wn_spc is the spacing, ie second coefficient
   # wco2 offset : delta_2 = 0.48 * ( delta_1 + wn_spc) - wn_spc)
   # sco2_offset : delta_3 = 0.38 * ( delta_1 + wn_spc) - wn_spc)

   disp_offsets = []
   
   for band_idx, band_scaling in enumerate(DISPERSION_OFFSET_SCALING):
      # Assuming 2 coefficients for describing dispersion for GOSAT
      # Since we need both to do modifcation to offset
      band_wn_start, band_spacing = dispersion_coefs[band_idx, :]

      # Use a different equation for l1b builds before v02.04.06
      if (l1b_build_id != None and (0,0,0) > l1b_build_id < (2,4,6)) or (index_scheme != None and index_scheme == 0):
         # For 0 based indexing in L1B file
         disp_obj.header['scaling_equation'] = '"0 based indexing"'
         curr_offset = DISPERSION_OFFSET_SCALING[band_idx] * (wn_shift - band_spacing) + band_spacing
      else:
         # For 1 based indexing in L1B file
         disp_obj.header['scaling_equation'] = '"1 based indexing"'
         curr_offset = DISPERSION_OFFSET_SCALING[band_idx] * wn_shift
      
      disp_offsets.append(curr_offset)

      # First column is index of parameter
      disp_obj.data[0, band_idx+1] = band_wn_start + curr_offset
      disp_obj.data[1, band_idx+1] = band_spacing

   disp_obj.header['offsets'] = disp_offsets
   disp_obj.header['solar_doppler_velocity'] = solar_doppler_velocity
   disp_obj.labels = DISPERSION_LABELS
   disp_obj.file_id = 'Dispersion apriori for: %s' % sounding_id
   disp_obj.write(apriori_out_file)

def create_dispersion_from_ascii(l1b_file, out_disp_file, disp_in_file, sounding_id=None, index_scheme=None):

    asc_l1b_obj = OCO_Matrix(l1b_file)

    if sounding_id == None:
        sounding_id = asc_l1b_obj.header['sounding_id']

    if disp_in_file == None:
        raise IOError('No dispersion file specified')

    latitude = float(asc_l1b_obj.header['sounding_latitude'].split()[0])
    sza_r    = math.radians(float(asc_l1b_obj.header['sounding_solar_zenith'].split()[0]))
    saz_r    = math.radians(float(asc_l1b_obj.header['sounding_solar_azimuth'].split()[0]))

    time_stamp = asc_l1b_obj.header['frame_time_stamp']
    time_struct = OCO_TextUtils.convert_timestamp_to_struct(time_stamp)

    pixel_ranges = asc_l1b_obj.pixel_ranges()
    aband_data = asc_l1b_obj['Radiance'][slice(*pixel_ranges[0]), 0]

    disp_in_obj = OCO_Matrix(disp_in_file)
    dispersion_coefs = disp_in_obj[DISPERSION_ASCII_COLUMN_IDENT].transpose()
    
    create_scene_dispersion_file(sounding_id, latitude, sza_r, saz_r, time_struct, aband_data, dispersion_coefs, out_disp_file, index_scheme=index_scheme)
    
def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
    if len(moduleSections) > 1:
        raise Exception('Only one instance of %s allowed per FILE block' % os.path.basename(__file__))

    l1b_file      = source
    out_disp_file = destination

    sounding_id  = Apply_Template(moduleSections[0].Get_Keyword_Value('sounding_id'), valuesDict, mapDict=mapDict)
    disp_in_file = Apply_Template(moduleSections[0].Get_Keyword_Value('dispersion_in_file'), valuesDict, mapDict=mapDict)
    index_scheme = Apply_Template(moduleSections[0].Get_Keyword_Value('index_scheme'), valuesDict, mapDict=mapDict)

    if index_scheme != None and len(index_scheme) > 0:
       index_scheme = int(index_scheme)

    if type(l1b_file) is str and h5py.is_hdf5(l1b_file):
        raise Exception('%s does not yet handle HDF L1B processing directly' % os.path.basename(__file__))
    else:
        create_dispersion_from_ascii(l1b_file, out_disp_file, disp_in_file, sounding_id, index_scheme)

def standalone_main():
   from optparse import OptionParser

   # Load command line options
   parser = OptionParser(usage="usage: %prog <spectra> <dispersion_in>...")

   parser.add_option( "-o", "--output_file", dest="output_file",
                      metavar="FILE",
                      default='./new_dispersion.dat',
                      help="where to write dispersion file other than default")

   parser.add_option( "-s", "--sounding_id", dest="sounding_id",
                      metavar="ID",
                      help="sounding id of dispersion information")

   parser.add_option( "-z", "--zero_indexing", dest="zero_indexing",
                      action="store_true",
                      help="force zero base indexing")

   # Parse command line arguments
   (options, args) = parser.parse_args()

   if (len(args) < 2):
      parser.error('Must supply spectra filename and dispersion input filename')

   l1b_file     = args[0]
   disp_in_file = args[1]

   if options.zero_indexing:
      index_scheme = 0
   else:
      index_scheme = 1

   create_dispersion_from_ascii(l1b_file, options.output_file, disp_in_file, options.sounding_id, index_scheme)

if __name__ == "__main__":
    standalone_main()
