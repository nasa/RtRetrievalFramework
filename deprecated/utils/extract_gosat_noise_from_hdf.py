#!/usr/bin/env python

import os
import math
import sys
import bisect
from optparse import OptionParser

import ACOS_File
import numpy

from OCO_Matrix import OCO_Matrix

COLOR_COLUMN_TMPL = 'COLOR_%d'
CNV_COLUMN_TMPL   = 'CNV_COEF_%d'

BAND_COLUMN_LABELS = []
for idx in range(len(ACOS_File.BAND_DATA_NAMES)):
    BAND_COLUMN_LABELS.append(COLOR_COLUMN_TMPL % (idx+1))
    BAND_COLUMN_LABELS.append(CNV_COLUMN_TMPL % (idx+1))

OUT_FILE_PREFIX = 'noise_cnv_%s_%s.dat'

NOISE_VALUES_HEADER_NAME = 'band_noise_values'
SOURCE_FILE_HEADER_NAME = 'source_acos_l1b_file'
SOURCE_SOUNDING_HEADER_NAME = 'source_sounding_id'

BAND_PARAM_HEADER_TMPL = 'band_modify_params_%d'

# Average top N radiances to derive signal
SIGNAL_AVG_COUNT = 10

FILE_ID = 'GOSAT Noise for Sounding: %s, Gain: %s.'

def extract_sounding_noise(l1b_obj, sounding_id, gain_code=None, modification_params=None):

    # Get the gain code used for current sounding id
    if gain_code == None:
        gain_code = l1b_obj.get_sounding_info('gain', sounding_id)
        
    out_obj = OCO_Matrix()

    if hasattr(gain_code, 'shape'):
        if numpy.all(gain_code == gain_code[0]):
            gain_code = gain_code[0]
        else:
            raise ValueError('Sounding %s does not have same gain code for both polarizations: %s' % (sounding_id, gain_code))
        
    noise_data = l1b_obj.get_error_data(sounding_id, calculate_noise=False, gain_code=gain_code)

    max_diodes = 0
    for band_cnv, band_noise_val in noise_data:
        num_diodes = band_cnv.shape[-1]
        max_diodes = max(max_diodes, num_diodes)

    if modification_params != None:
        radiances = l1b_obj.get_radiance_data(sounding_id)

    band_noise_values = []
    cnv_data = numpy.zeros((max_diodes, len(BAND_COLUMN_LABELS)), dtype=float)
    for band_idx, band_noise_tuple in enumerate(noise_data):
        band_cnv, band_noise_val = band_noise_tuple

        if len(band_cnv.shape) > 1:
            band_cnv = numpy.average( band_cnv, axis=0 )

        if hasattr(band_noise_val, 'shape') and band_noise_val.shape[0] > 1:
            band_noise_val = numpy.average( band_noise_val, axis=0 )
        
        num_diodes = band_cnv.shape[0]
        
        band_noise_values.append( band_noise_val)

        if modification_params != None:
            band_params = modification_params[band_idx]
            out_obj.header[BAND_PARAM_HEADER_TMPL % (band_idx+1)] = ' '.join([str(curr_val) for curr_val in band_params])

            # Derive signal from average of top N radiance values
            band_radiance = radiances[band_idx]

            if len(band_radiance.shape) > 1:
                band_radiance = numpy.average(band_radiance, axis=0)
            
            band_radiance.sort()
            signal = numpy.mean(band_radiance[-SIGNAL_AVG_COUNT:])

            # Calculated updated band_cnv
            param_a1, param_b1 = band_params
            band_cnv[:] = (param_a1 * signal) / band_noise_val + param_b1 * band_cnv[:]
            
        cnv_data[:num_diodes, band_idx*2]   = numpy.arange(1,num_diodes+1)
        cnv_data[:num_diodes, band_idx*2+1] = band_cnv[:]

    out_obj.header[NOISE_VALUES_HEADER_NAME] = ' '.join([ ('%e' % val) for val in band_noise_values ])
    out_obj.header[SOURCE_FILE_HEADER_NAME] = os.path.realpath(l1b_obj.filename)
    out_obj.header[SOURCE_SOUNDING_HEADER_NAME] = sounding_id
    out_obj.file_id = FILE_ID % (sounding_id, gain_code)
    out_obj.labels = BAND_COLUMN_LABELS
    out_obj.data = cnv_data
    
    return out_obj

def extract_all_noise(hdf_filename, output_directory, sounding_id=None, modify_file=None):

    modification_params=None
    if modify_file != None:
        modify_obj = OCO_Matrix(modify_file)
        modification_params = []
        for col_idx in range(modify_obj.dims[1]):
            modification_params.append( modify_obj.data[:,col_idx] )

    with ACOS_File.L1B(hdf_filename) as l1b_obj:
        if sounding_id == None:
            sounding_id = l1b_obj.get_sounding_ids()[0]
            print >>sys.stderr, 'Using just the first sounding id: %s' % sounding_id

        for gain_code, gain_tmpl in ACOS_File.GOSAT_CNV_COEF_DATASET.items():
            output_filename = os.path.join(output_directory, OUT_FILE_PREFIX % (sounding_id, gain_code))
            file_obj = extract_sounding_noise(l1b_obj, sounding_id, gain_code, modification_params)
            print 'Writing: %s' % output_filename
            file_obj.write(output_filename)

def standalone_main():
    # Load command line options
    parser = OptionParser(usage="usage: %prog [options] <hdf_filename>")

    parser.add_option( "-o", "--output_dir", dest="output_dir",
                       metavar="DIR", default="./",
                       help="directory where extracted files are written")

    parser.add_option( "-s", "--sounding_id", dest="sounding_id",
                       metavar="STR", 
                       help="desired sounding id for extraction")

    parser.add_option( "-m", "--modify_file", dest="modify_file",
                       metavar="STR", 
                       help="empirical noise correction file")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('Must supply hdf filename')

    extract_all_noise(args[0], options.output_dir, options.sounding_id, options.modify_file)


if __name__ == "__main__":
    standalone_main()
