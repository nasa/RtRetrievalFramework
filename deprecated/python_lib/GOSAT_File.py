# Read GOSAT TANSO FTS L1B file and add a programatic wrapper interface

import os
import sys
import time
import datetime

import numpy
import h5py

import OCO_MathUtil
from OCO_Matrix import OCO_Matrix

SWIR_BANDS  = (1,2,3)
POL_ORDER   = ('P', 'S')

# List of tuples where first index is 1d index for packed BAND*POL datasets
# second index is band index in 2d BANDxPOL matrix
# third index is pol index in 2d BANDxPOL matrix"""

PACKED_BAND_POL_INDEXES = [ (x*len(POL_ORDER)+y, x,y) for x in range(len(SWIR_BANDS)) for y in range(len(POL_ORDER)) ]

EXPOSURE_START_TIME_DATASET  = '/exposureAttribute/pointAttribute/Time'
ZPD_TIME_DATASET             = '/exposureAttribute/pointAttribute/RadiometricCorrectionInfo/ZPDPassTime_SWIR'

OBS_LEN_DATASET = '/exposureAttribute/numPoints'

SWIR_SPEC_DATASET_TMPL = 'Spectrum/SWIR/band%d/obsWavelength'
SWIR_DISP_DATASET      = '/exposureAttribute/pointAttribute/RadiometricCorrectionInfo/spectrumObsWavelengthRange_SWIR'

DISP_CORRECTIONS = ( -0.4326, -0.2078, -0.1621 )

RAD_CNV_FILE_TMPL     = os.path.join( os.environ['L2_INPUT_PATH'], 'gosat/radiometric/TANSO-FTS-RcalXX0620_B%d.dat')
RAD_CNV_UNIT_COL      = 'Wavenumber'
RAD_CNV_SENS_COL_TMPL = 'CNV(B%d%s_gain%s)'

GAIN_DATASET = '/exposureAttribute/pointAttribute/RadiometricCorrectionInfo/gain_SWIR'

GEOMETRY_DATASETS = {
    'latitude': '/exposureAttribute/pointAttribute/geometricInfo/centerLat',
    'longitude': '/exposureAttribute/pointAttribute/geometricInfo/centerLon',
    'azimuth':   '/exposureAttribute/pointAttribute/geometricInfo/AZ',
    'zenith':    '/exposureAttribute/pointAttribute/geometricInfo/AOI',
    'solar_azimuth': '/exposureAttribute/pointAttribute/sun/AZ',
    'solar_zenith': '/exposureAttribute/pointAttribute/sun/EL',
    }

class L1B(h5py.File):

    def __init__(self, filename, sensitivity_correction=True, correct_disp=False):
        h5py.File.__init__(self, filename, 'r')

        self.sens_corr_data = None
        if sensitivity_correction:
            try:
                self._load_sens_conv_files()
            except IOError as err:
                print >>sys.stderr, 'Could not load sensitivity correction data, will not make correction: %s' % err

        self.correct_disp = correct_disp

    def get_num_observations(self):
        return self[OBS_LEN_DATASET][0]

    def get_exposure_time(self, obs_index):
        self._check_obs_index(obs_index)
        
        return self[EXPOSURE_START_TIME_DATASET][obs_index]
   
    def get_observation_id(self, obs_index):

        time_values = self.get_exposure_time(obs_index)

        id_format = list(time_values[0:5]) + [round(time_values[5])]
        obs_id = '%04d%02d%02d%02d%02d%02d' % tuple(id_format)

        return numpy.array(numpy.int64(obs_id))

    def get_time_stamp(self, obs_index):

        time_stamp = '%04d-%02d-%02dT%02d:%02d:%06.3fZ' % self.get_exposure_time(obs_index)

        return time_stamp

    def get_epoch_times(self, obs_index):
        exposure_start = self.get_exposure_time(obs_index)

        epoch_times = numpy.zeros((len(SWIR_BANDS), len(POL_ORDER)), dtype=float)
        for packed_idx, band_idx, pol_idx in PACKED_BAND_POL_INDEXES:
            zpd_tuple = list(exposure_start)

            zpd_seconds = self[ZPD_TIME_DATASET][obs_index, packed_idx]

            dt_inp = list(zpd_tuple[:5]) + [ int(zpd_tuple[5]), int((zpd_tuple[5] - int(zpd_tuple[5]))*1e6) ]
            epoch_times[band_idx, pol_idx] = time.mktime((datetime.datetime(*dt_inp) + datetime.timedelta(seconds=zpd_seconds)).utctimetuple())

        return epoch_times
    def get_oco_times(self, obs_index):
        zpd_times = self.get_epoch_times(obs_index)

        oco_start_epoch = time.mktime(datetime.datetime(1993,1,1).utctimetuple())

        for band_idx in range(len(SWIR_BANDS)):
            for pol_idx in range(len(POL_ORDER)):
                zpd_times[band_idx][pol_idx] -= oco_start_epoch
                            
        return zpd_times

    def get_swir_geometry(self, obs_index):

        geom_dict = {}
        for dest_name, dataset in GEOMETRY_DATASETS.items():
            if len(self[dataset].shape) == 1:
                geom_dict[dest_name] = numpy.array(self[dataset][obs_index])
            else:
                geom_matrix = numpy.zeros((len(SWIR_BANDS), len(POL_ORDER)), dtype=float)
                for packed_idx, band_idx, pol_idx in PACKED_BAND_POL_INDEXES:
                    geom_matrix[band_idx, pol_idx] = self[dataset][packed_idx, obs_index]
                geom_dict[dest_name] = geom_matrix

        # Adjust angles
        if geom_dict.has_key('azimuth'):
            geom_dict['azimuth'] = 180.0 - geom_dict['azimuth']

        if geom_dict.has_key('zenith'):
            geom_dict['zenith'] = 90.0 + geom_dict['zenith']

        return geom_dict

    def get_swir_spectrum(self, obs_index):

        self._check_obs_index(obs_index)
                
        swir_spectrum = []
        for curr_band in SWIR_BANDS:
            # Get obs real S and P spectrum
            band_wns = self.get_swir_wavenumbers(obs_index)[SWIR_BANDS.index(curr_band)]
            band_data = numpy.transpose(self[SWIR_SPEC_DATASET_TMPL % curr_band][obs_index, :, :, 0])
            gain = [chr(gain_int) for gain_int in self[GAIN_DATASET][obs_index, :, 0]]
            
            band_data = self._apply_sens_correction(curr_band, band_wns, band_data, gain)
            swir_spectrum.append(band_data)

        return swir_spectrum

    def get_swir_num_wavenumbers(self):
        num_wn = []
        for curr_band in SWIR_BANDS:
            num_wn.append( self[SWIR_SPEC_DATASET_TMPL % curr_band].shape[2] )

        return num_wn

    def get_swir_wavenumbers(self, obs_index):

        self._check_obs_index(obs_index)

        num_wn = self.get_swir_num_wavenumbers()
            
        wavenumbers = []
        for band_index, dataset_band_index in enumerate(range(0,len(SWIR_BANDS)*2, 2)):
            band_coefs_p1 = self[SWIR_DISP_DATASET][obs_index, dataset_band_index, :]
            band_coefs_p2 = self[SWIR_DISP_DATASET][obs_index, dataset_band_index+1, :]

            if self.correct_disp:
                band_coefs_p2 += DISP_CORRECTIONS[band_index]

            band_wn = numpy.zeros((num_wn[band_index], 2), dtype=float)
            band_wn[:, 0] = band_coefs_p1[1] + band_coefs_p1[0] * numpy.arange(1,num_wn[band_index]+1)
            band_wn[:, 1] = band_coefs_p2[1] + band_coefs_p2[0] * numpy.arange(1,num_wn[band_index]+1)
            wavenumbers.append(band_wn)

        return wavenumbers

    def _check_obs_index(self, obs_index):
        if obs_index < 0 or obs_index >= self.get_num_observations():
            raise Exception('observation index supplied: %d out of range [0, %d]' % (obs_index, self.get_num_observations()))

    def _load_sens_conv_files(self):

        self.sens_corr_data = []
        for curr_band in SWIR_BANDS:
            corr_file = RAD_CNV_FILE_TMPL % curr_band
            corr_obj  = OCO_Matrix(corr_file, ignore_conv_err=True)
            self.sens_corr_data.append(corr_obj)

    def _apply_sens_correction(self, band_num, band_wns, band_data, gain):

        if self.sens_corr_data == None:
            return band_data

        band_index = SWIR_BANDS.index(band_num)

        corr_obj = self.sens_corr_data[band_index]

        unit_col_idx = corr_obj.labels.index(RAD_CNV_UNIT_COL)

        for pol_idx, pol_name in enumerate(POL_ORDER):
            corr_col_name = RAD_CNV_SENS_COL_TMPL % (band_num, pol_name, gain[pol_idx])

            try:
                corr_col_idx = corr_obj.labels.index(corr_col_name)
            except:
                raise IOError('Could not find column named: %s in sensitivity file: %s' % (corr_col_name, corr_obj.filename))


            # Retrieve correction data
            resampled_corr_data = numpy.zeros(band_data.shape[0], dtype=float)

            src_corr_data = corr_obj.data[:,corr_col_idx]
            resampled_corr_data = OCO_MathUtil.linear_interpol(src_corr_data, corr_obj.data[:, unit_col_idx], band_wns[:,pol_idx], extrapolate=False, extrap_err=False)
            
            valid_points = numpy.where(abs(resampled_corr_data) > 0.0)[0]
            band_data[valid_points, pol_idx] *= resampled_corr_data[valid_points]

        return band_data
