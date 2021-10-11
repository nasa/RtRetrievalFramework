from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from past.utils import old_div
import six
# Read ACOS,OCO L1B file and add a programatic wrapper interface

import os
import re
import sys
import time
import copy
import types
import bisect
import datetime
import traceback
from operator import itemgetter
from collections import namedtuple, Counter
import logging

import numpy
import h5py
from . import tai64n
import six

logger = logging.getLogger()

# Require at least h5py version 1.2.x
if h5py.version.version_tuple < (1,2,0):
    raise ImportError('At least version 1.2.0 of h5py is required, you have: %s' % h5py.version.version)


from .surface_prop import ModisEcoMap

GOSAT_INST_NAME = 'gosat'
OCO_INST_NAME   = 'oco'

SOUNDING_ID_DATASET = { OCO_INST_NAME:   'SoundingGeometry/sounding_id',
                        GOSAT_INST_NAME: 'SoundingHeader/sounding_id'
                        }

SOUNDING_ID_DIMENSIONS = { SOUNDING_ID_DATASET[OCO_INST_NAME]:   ('Frame', 'Sounding'),
                           SOUNDING_ID_DATASET[GOSAT_INST_NAME]: ('Exposure', 'Polarization'), }

RADIANCE_GROUP    = { OCO_INST_NAME:   'SoundingMeasurements',
                      GOSAT_INST_NAME: 'SoundingSpectra' }

BAND_DATA_NAMES   = ('o2', 'weak_co2', 'strong_co2')
RADIANCE_DATASETS = [ 'radiance_%s' % data_name for data_name in BAND_DATA_NAMES]

MAX_MEAS_SIGNAL = { OCO_INST_NAME:   (1.4e21, 4.9e20, 1.7e20),
                    GOSAT_INST_NAME: (7.2e-6, 6.5e-6, 4.5e-6),
                    }

OCO_SNR_COEF_DATASET = '/InstrumentHeader/snr_coef'

GOSAT_CNV_COEF_DATASET = { 'H': '/InstrumentHeader/cnv_coef_highgain_%s',
                           'M': '/InstrumentHeader/cnv_coef_medgain_%s',
                           }
GOSAT_NOISE_DATASET = '/SoundingSpectra/noise_%s'

TIME_STRING_FORMAT = '%Y-%m-%dT%H:%M:%S'
TIME_STRING_EXPECT_LEN = len(time.strftime(TIME_STRING_FORMAT))
TIME_DATASET = ('/FootprintGeometry/footprint_time',
                '/FootprintGeometry/footprint_time_string',
                '/FootprintGeometry/footprint_time_tai93',
                )

# These names match what is needed to create a FP sounding info file for
# simulations or that need to be present in an ASCII L1B file
SOUNDING_INFO_DATASETS = { 'frame_time_stamp':       { OCO_INST_NAME:   '/FrameHeader/frame_time_string',
                                                       GOSAT_INST_NAME: '/SoundingHeader/exposure_start_time_string' },
                           'polarization_angle':     { OCO_INST_NAME:   '/FootprintGeometry/footprint_polarization_angle',
                                                       GOSAT_INST_NAME: None,
                                                       },
                           'stokes_coefficients':    { OCO_INST_NAME:   'FootprintGeometry/footprint_stokes_coefficients',
                                                       GOSAT_INST_NAME: ('/FootprintGeometry/footprint_mueller_matrix',
                                                                         '/FootprintGeometry/stokes_coefficients',
                                                                         '/FootprintGeometry/footprint_stokes_coefficients',
                                                                         ),
                                                       },
                           'relative_velocity':      { OCO_INST_NAME:   '/FrameGeometry/relative_velocity',
                                                       GOSAT_INST_NAME: '/SpacecraftGeometry/relative_velocity',
                                                       },
                           'sounding_latitude':      '/FootprintGeometry/footprint_latitude',
                           'sounding_longitude':     '/FootprintGeometry/footprint_longitude',
                           'sounding_altitude':      '/FootprintGeometry/footprint_altitude',
                           'sounding_azimuth':       '/FootprintGeometry/footprint_azimuth',
                           'sounding_zenith':        '/FootprintGeometry/footprint_zenith',
                           'sounding_solar_azimuth': '/FootprintGeometry/footprint_solar_azimuth',
                           'sounding_solar_zenith':  '/FootprintGeometry/footprint_solar_zenith',
                           'spacecraft_lat':         { OCO_INST_NAME:   '/FrameGeometry/spacecraft_lat',
                                                       GOSAT_INST_NAME: '/SpacecraftGeometry/spacecraft_lat',
                                                       },
                           'spacecraft_lon':         { OCO_INST_NAME:   '/FrameGeometry/spacecraft_lon',
                                                       GOSAT_INST_NAME: '/SpacecraftGeometry/spacecraft_lon',
                                                       },
                           'spacecraft_alt':         { OCO_INST_NAME:   '/FrameGeometry/spacecraft_alt',
                                                       GOSAT_INST_NAME: '/SpacecraftGeometry/spacecraft_alt',
                                                       },
                           'land_fraction':          { OCO_INST_NAME:   None,
                                                       GOSAT_INST_NAME: '/SoundingGeometry/sounding_land_fraction',
                                                       },

                           'land_water_indicator':   { OCO_INST_NAME:   '/SoundingGeometry/sounding_land_water_indicator',
                                                       GOSAT_INST_NAME: None,
                                                       },
                           }
# Additional info datasets that wount be put into a header for an L1B ascii file
ADDL_INFO_DATASETS    = { 'gain':       { OCO_INST_NAME:    None,
                                          GOSAT_INST_NAME:  ('SoundingHeader/gain_SWIR', 'SoundingHeader/gain_swir') },
                          'dispersion': { OCO_INST_NAME:    ('/Metadata/DispersionCoefSamp', '/InstrumentHeader/dispersion_coef_samp'),
                                          GOSAT_INST_NAME:  '/SoundingHeader/wavenumber_coefficients',
                                          },
                          }

# Set up alternative names for quicker access
for snd_name, dataset in list(SOUNDING_INFO_DATASETS.items()):
    if snd_name.find('sounding_') == 0:
        ADDL_INFO_DATASETS[snd_name.replace('sounding_', '')] = SOUNDING_INFO_DATASETS[snd_name]

#########
# All of the following set up default and specific shape names when none available
DEFAULT_L1B_INFO_SHAPE_NAMES = { OCO_INST_NAME: ('Frame', 'Sounding', 'Band'),
                                 GOSAT_INST_NAME: ('Exposure', 'Band', 'Polarization'),
                                 }

# Set up specific dataset shapes
SPECIFIC_L1B_INFO_SHAPE_NAMES = { ADDL_INFO_DATASETS['dispersion'][GOSAT_INST_NAME]: ('Exposure','Band','Polarization','WNConv'),
                                  ADDL_INFO_DATASETS['dispersion'][OCO_INST_NAME]: ('Spectrum','Sounding', 'DispersionCoefficient'),
                                  SOUNDING_INFO_DATASETS['frame_time_stamp'][GOSAT_INST_NAME]: ('Exposure',),
                                  SOUNDING_INFO_DATASETS['frame_time_stamp'][OCO_INST_NAME]: ('Frame',),
                                  SOUNDING_INFO_DATASETS['relative_velocity'][GOSAT_INST_NAME]: ('Exposure',),
                                  SOUNDING_INFO_DATASETS['relative_velocity'][OCO_INST_NAME]: ('Frame',),
                                  SOUNDING_INFO_DATASETS['land_fraction'][GOSAT_INST_NAME]: ('Exposure',),
                                  SOUNDING_INFO_DATASETS['land_water_indicator'][OCO_INST_NAME]: ('Frame', 'Sounding'),

                                  }
# Set shapes for gain
for gain_ds in ADDL_INFO_DATASETS['gain'][GOSAT_INST_NAME]:
    SPECIFIC_L1B_INFO_SHAPE_NAMES[gain_ds] = ('Exposure', 'Polarization')

# Set shapes of named band datasets
for band_name in BAND_DATA_NAMES:
    SPECIFIC_L1B_INFO_SHAPE_NAMES['/%s/noise_%s' % (RADIANCE_GROUP[GOSAT_INST_NAME], band_name)] = ('Exposure', 'Polarization')
    SPECIFIC_L1B_INFO_SHAPE_NAMES['/%s/radiance_%s' % (RADIANCE_GROUP[GOSAT_INST_NAME], band_name)] = ('Exposure', 'Polarization', 'Color')
    
    for gain_code, conv_ds_tmpl in list(GOSAT_CNV_COEF_DATASET.items()):
        SPECIFIC_L1B_INFO_SHAPE_NAMES[GOSAT_CNV_COEF_DATASET[gain_code] % band_name] = SPECIFIC_L1B_INFO_SHAPE_NAMES['/SoundingSpectra/radiance_%s' % band_name]

# Set shapes for spacecraft datasets
for shape_name, ds_desc in list(SOUNDING_INFO_DATASETS.items()):
    if shape_name.find('spacecraft_') == 0:
        SPECIFIC_L1B_INFO_SHAPE_NAMES[SOUNDING_INFO_DATASETS[shape_name][GOSAT_INST_NAME]] = ('Exposure',)
        SPECIFIC_L1B_INFO_SHAPE_NAMES[SOUNDING_INFO_DATASETS[shape_name][OCO_INST_NAME]] = ('Frame',)


# Set stokes shape names
for stokes_ds in SOUNDING_INFO_DATASETS['stokes_coefficients'][GOSAT_INST_NAME]:
    SPECIFIC_L1B_INFO_SHAPE_NAMES[stokes_ds] = ('Exposure', 'Band', 'Polarization', 'StokesCoefficient')

#########
# Set build dataset name
BUILD_ID_DATASET = { OCO_INST_NAME:   '/Metadata/BuildId',
                     GOSAT_INST_NAME: '/Metadata/BuildId',
                     }

GOSAT_POL_ORDER   = ('P', 'S')

# When TAI93 time starts
TAI93_START  = datetime.datetime(1993,1,1,0,0)
TAI93_LPSECS = tai64n.__tai_seconds(TAI93_START)

LAND_FRACTION_PERCENTAGE = { # Consistent with SDOS sounding picking procedures
                             'land':  lambda lf: 90.0 <  lf <= 100.0,
                             
                             # for now classify these as ocean so populator
                             # does not choke when someone feeds it sounding
                             # list that is not all land
                             'ocean': lambda lf: 0.0  <= lf <=  90.0,
                             }

# 0 -> Land
# 1 -> Ocean
# 2 -> Inland water
# 3 -> Mixed land water
LAND_WATER_INDICATOR = { 'land': lambda li: li == 0 or li == 3,
                         'ocean': lambda li: li == 1 or li == 2,
                         }
########

class AveragingError(ValueError):
    "Trouble peforming averaging of data from L1B file"

class MissingDataset(LookupError):
    "Dataset not found in HDF file"

class MissingSounding(LookupError):
    "Sounding not found in HDF file"

########

class NamedShapeArray(numpy.ndarray):
    """Stores a a view of a numpy array where the shapes of that data
    are described by a named tuple."""
        
    def __new__(cls, data, indexes, preserve_shape=False, **kwargs):
        if not isinstance(indexes, tuple) and hasattr(indexes, "_fields"):
            raise TypeError("expected a named tuple describing shape names and indexes to extract from data")

        new_dim_names = []
        extract_indexes = []
        for dim_name, dim_extract in zip(indexes._fields, tuple(indexes)):            
            if isinstance(dim_extract, int):
                if preserve_shape:
                    # Preserve shape by converting any integer indexing to a slice
                    new_dim_names.append(dim_name)
                    extract_indexes.append(slice(dim_extract,dim_extract+1))
                else:
                    # Ignore this dim name since it will be reduced
                    extract_indexes.append(dim_extract)
            else:
                new_dim_names.append(dim_name)
                extract_indexes.append(dim_extract)

        # Do things this way so that indexing of the result
        # turns out as expected
        try:
            if len(extract_indexes) == 1 and not isinstance(extract_indexes[0], slice):
                buf_data = data[numpy.array(extract_indexes[0])]
            else:
                buf_data = data[tuple(extract_indexes)]
        except MemoryError as e:
            raise MemoryError("Ran out of memory initializing NamedShapeArray")

        new_obj = numpy.asarray(buf_data).view(cls)

        # Make sure we have an array and not something reduced to scalar
        new_obj.inp_dim_names = new_dim_names
        data_shp_class = namedtuple('shape', new_dim_names)
        if hasattr(buf_data, 'shape') and len(new_dim_names) == len(buf_data.shape):
            new_obj.named_shape = data_shp_class(*buf_data.shape)
        elif len(new_dim_names) == 0:
            new_obj.named_shape = data_shp_class()
        else:
            new_obj.named_shape = None

        return new_obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.named_shape = getattr(obj, 'named_shape', None)

def is_list_like(d):
    '''Check if an object is like a list. Basically we check for __iter__,
    but explicitly exclude str types which in python 2 did not have __iter__,
    but was added in python 3 (our code use to just look for __iter__, but
    this breaks in python 3 since both ["foo"] and "foo" have __iter__.'''
    if isinstance(d, six.string_types):
        return False
    return hasattr(d, '__iter__')

################################################################################################

class SoundingDataFile(h5py.File):
    def __init__(self, filename, mode, **kwargs):
        # Init super-class
        try:
            if mode == None:
                mode = 'r'
            h5py.File.__init__(self, filename, mode, **kwargs)
        except IOError as exc:
            raise IOError("Could not open hdf file: %s with mode: %s, due to error: %s" % (filename, mode, exc))

        self._sounding_id_dataset = None
        self._data_shape_name_dict = {}
        self._default_shape_names = ()

        # Cache these so a re-read of the file is unnecessary 
        self._id_dim_names = None

    def get_sounding_ids(self):
        try:
            file_sounding_ids = self[self._sounding_id_dataset]
        except KeyError:
            raise KeyError("Unable to find sounding id dataset: %s in file: %s" % (self._sounding_id_dataset, self.filename))
        return file_sounding_ids

    def get_id_dim_names(self):
        if self._id_dim_names == None:
            self._id_dim_names = self.get_data_shape(self._sounding_id_dataset)
        return self._id_dim_names 

    def get_sounding_indexes(self, sounding_id):
        """Find sounding id through bisection, assumes dataset is flat"""

        # If not already present get a flattened list of ids
        # and their n-dimensional index
        if not hasattr(self, "_sounding_id_to_index"):
            sounding_ids = self.get_sounding_ids()
            self._id_type = sounding_ids.dtype.type
            
            self._sounding_id_to_index = dict((sid, index) \
                for (index, sid) in numpy.ndenumerate(sounding_ids))

        # Create a named tuple using the names of the dimensions for the indexes
        dim_names = self.get_id_dim_names()
        index_class = namedtuple('indexes', dim_names)

        # Find the sounding id in the flattened dataset and get its n dim index
        id_indexes = self._sounding_id_to_index[self._id_type(sounding_id)]

        # Determine if we actually found the value
        if id_indexes is None:
            raise MissingSounding('Did not find correct sounding id: %s in hdf file: %s' % (sounding_id, self.filename))

        # Do some sanity checking
        try:
            found_snd_id = self.get_sounding_ids()[id_indexes]
        except ValueError:
            raise ValueError('Indexes: %s for sounding id: %s out of range for sounding id matrix of shape: %s' % (id_indexes, sounding_id, sounding_ids.shape))
        
        # Construct a named tuple with those indexes, padding None if
        # the named tuple has more names than indexes we found
        if not is_list_like(id_indexes):
            id_indexes = [ id_indexes ]
        else:
            id_indexes = list(id_indexes)

        while len(id_indexes) < len(dim_names):
            id_indexes.append(None)

        index_tuple = index_class(*id_indexes[:len(dim_names)])

        return index_tuple

    def get_data_shape(self, dataset, dflt_dims_map=SOUNDING_ID_DIMENSIONS):
        """Try to determine shape names either from HDF file itself or using
        default names if none were passed explicitly"""
       
        if isinstance(dataset, h5py.Dataset):
            dataset_obj = dataset
        else:
            dataset_obj = self[dataset]

        dataset_shape = dataset_obj.attrs.get('Shape', None)

        if dataset_shape is not None:
            # Intelligently split up the shape name so that something like Dim_1 does
            # not get split incorrectly 
            # Ensure we are using a unicode string instead of bytes as HDF lib will return
            shape_names = []

            # Work arond differences between h5py versions
            try:
                dataset_shape = dataset_shape[0].replace("_Array", '')
            except:
                dataset_shape = dataset_shape[0].decode('utf-8').replace("_Array", '')

            # Search for next location of _
            while dataset_shape.find("_") >= 0:
                split_loc = dataset_shape.find("_")
                end_loc = split_loc
                check_num_idx = split_loc + 2 
                # While the portion after _ is still a number keep increasing the end location
                while dataset_shape[split_loc+1:check_num_idx].isdigit() and check_num_idx <= len(dataset_shape):
                    check_num_idx += 1
                    end_loc = check_num_idx - 1

                # Record new shape name and prune it off shape string
                shape_names.append(dataset_shape[:end_loc])
                dataset_shape = dataset_shape[end_loc+1:]

            # Append any final string
            if len(dataset_shape) > 0:
                shape_names.append(dataset_shape)

        elif len(self._data_shape_name_dict) > 0 and dataset_obj.name in self._data_shape_name_dict:
            shape_names = self._data_shape_name_dict[dataset_obj.name]
        elif len(self._default_shape_names) == len(dataset_obj.shape):
            # Use only if shapes match
            shape_names = self._default_shape_names
        elif dflt_dims_map != None and dataset_obj.name.strip('/') in dflt_dims_map:
            shape_names = dflt_dims_map[dataset_obj.name.strip('/')]
        elif hasattr(dataset_obj, "shape"):
            # Create generic names and try and assign the id names where we can
            shape_names = ['Dim%d' % dim for dim in range(len(dataset_obj.shape))]

            snd_ids_shape = self.get_sounding_ids().shape
            id_shape_name = self.get_id_dim_names() 

            # Match the size of the datasets dimenisons to those in the sounding ids
            # Assign the appropriate id name when matching
            can_match_idx = 0
            for shape_idx, curr_name in enumerate(shape_names):
                dim_size = dataset_obj.shape[shape_idx] 
                if dim_size in snd_ids_shape:
                    which_id = snd_ids_shape.index(dim_size)
                    # Make sure we match the id names in order,
                    # to avoid mismatches with dimensions of the sam
                    # size as the smaller of the id dimensions
                    if which_id == can_match_idx and which_id < len(id_shape_name):
                        shape_names[shape_idx] = id_shape_name[which_id]
                        can_match_idx += 1
        else:
            shape_names = ()

        if len(shape_names) > 0 and len(shape_names) != len(dataset_obj.shape):
            raise ValueError("Size of dataset shape names: %s does not match size of data shape: %s" % (shape_names, dataset_obj.shape))

        return shape_names

    def get_sounding_data(self, data_name, sounding_id=None, indexes=None, return_indexes=False, **kwargs):
        """Retrieves data from the file given a partial name that is matched against
        the dataset names present in the file. A full dataset name can use either /
        for group seperator or __ for cases where you want the data_name to be a
        valid identifier.

        If sounding_id is supplied then only this sounding id's data is returned.

        If indexes is supplied then those indexes are used along the sounding dataset dimension.

        Both keyword arguments can not be supplied simultaneously."""

        # Replace any __ with / so inputted name can be a valid identifier
        # and still match the full dataset name
        data_name = data_name.replace("__", "/")

        # Try to see if data name is a full dataset name
        dataset_obj = self.select_valid_dataset(data_name, return_obj=True)

        # Match data_name against data sets, may match multiple
        if dataset_obj == None:
            dataset_obj = []
            for group_name, group_obj in self.items():
                for ds_name, ds_obj in group_obj.items():
                    dataset_name = "%s/%s" % (group_name, ds_name)
                    if re.search(data_name, dataset_name):
                        dataset_obj.append( ds_obj )

        # See if we matched too many or none
        if dataset_obj == None or is_list_like(dataset_obj) and len(dataset_obj) == 0:
            raise ValueError("No datasets matched data name: %s" % data_name)
        elif (is_list_like(dataset_obj) and not hasattr(dataset_obj, "shape")) and len(dataset_obj) > 1:
            raise ValueError("Data name: %s matches too many datasets: %s" % (data_name, [ o for o in dataset_obj]))
        elif is_list_like(dataset_obj) and not hasattr(dataset_obj, "shape"):
            dataset_obj = dataset_obj[0]

        data_dim_indexes = self.get_data_indexes(dataset_obj, sounding_id, indexes)

        try:
            named_shape = NamedShapeArray(dataset_obj[:], data_dim_indexes, **kwargs)
        except AttributeError as e:
            raise AttributeError(str(e) + " for dataset: %s" % dataset_obj.name)
        except ValueError as e:
            raise ValueError(str(e) + " for dataset: %s" % dataset_obj.name)
        
        return named_shape

    def get_data_indexes(self, dataset_obj, sounding_id=None, indexes=None):
        if sounding_id != None and indexes != None:
            raise Exception("Can not supply sounding id and indexes simultaneously")
        elif sounding_id != None and indexes == None:
            indexes = [self.get_sounding_indexes(sounding_id)]
        elif indexes == None:
            indexes = numpy.arange(len(self.get_sounding_ids()))

        shape_names = list(self.get_data_shape(dataset_obj.name))
        if shape_names == None:
            raise ValueError('No shape names are defined for dataset: %s in file: %s' % (dataset_name, self.filename))
        elif len(shape_names) != len(dataset_obj.shape):
            raise ValueError('Length of shape names: %s does not match length of data: %s for: %s' % (len(shape_names), len(dataset_obj.shape), dataset_obj.name))

        # Data shape name
        id_dims = self.get_id_dim_names()

        if len(id_dims) > 1:
            index_get = [ itemgetter(idx) for idx in range(len(id_dims)) ]
        else:
            index_get = [ lambda f:f ]

        data_dim_indexes = []
        for dim_idx, shp_name in enumerate(shape_names):
            if shp_name in id_dims:
                shape_idxs = tuple(map(index_get[id_dims.index(shp_name)], indexes))

                if len(shape_idxs) == 1:
                    shape_idxs = shape_idxs[0]
                data_dim_indexes.append(shape_idxs)
            else:
                data_dim_indexes.append( slice(dataset_obj.shape[dim_idx]) )

        # Eliminate duplicate shape names
        for name, num_occurance in list(Counter(shape_names).items()):
            if num_occurance > 1:
                for name_count in range(1, num_occurance+1):
                    shape_names[shape_names.index(name)] = "%s_%d" % (name, name_count)

        data_idx_class = namedtuple('indexes', tuple(shape_names))
        named_indexes = data_idx_class(*data_dim_indexes)
       
        return named_indexes

    def select_valid_dataset(self, dataset_names, return_obj=False):
        """Queries the HDF file for the first of multiple dataset names that are present within the file.
        Can optionally return the dataset object instead of the found dataset name.
        """

        # A iterable means we should try multiple names
        if is_list_like(dataset_names):
            for curr_dataset_try in dataset_names:
                dataset_obj = self.get(curr_dataset_try, None)
                found_dataset_name = curr_dataset_try

                if dataset_obj != None:
                    break
        else:
            found_dataset_name = dataset_names
            if return_obj:
                dataset_obj = self.get(found_dataset_name, None)

        if return_obj:
            return dataset_obj
        else:
            return found_dataset_name


################################################################################################

class L1B(SoundingDataFile):
    
    def __init__(self, filename, mode=None, **kwargs):
        # Init super-class
        SoundingDataFile.__init__(self, filename, mode, **kwargs)

        # Initialize as needed
        self.surf_type_obj = None

        # Needed internally for selecting correct dataset
        self.instrument_name = self.determine_instrument_name()

        # See if the file type is one we recognize
        if self.instrument_name not in (GOSAT_INST_NAME, OCO_INST_NAME):
            raise LookupError('Unrecognized instrument name detected: %s' % self.instrument_name)

        self._sounding_id_dataset = SOUNDING_ID_DATASET[self.instrument_name]

        self._data_shape_name_dict = SPECIFIC_L1B_INFO_SHAPE_NAMES
        self._default_shape_names = DEFAULT_L1B_INFO_SHAPE_NAMES[self.instrument_name]
        
    def get_sounding_ids(self, add_polarization=False):
        file_sounding_ids = SoundingDataFile.get_sounding_ids(self)
        
        if add_polarization and self.instrument_name == GOSAT_INST_NAME:
            sounding_ids = []

            for curr_id in file_sounding_ids:
                for pol_name in GOSAT_POL_ORDER:
                    sounding_ids.append( '%d%s' % (curr_id, pol_name) )

            return sounding_ids
        else:
            return file_sounding_ids

    def get_id_dim_names(self):
        dflt_dim_names = SOUNDING_ID_DIMENSIONS.get(self._sounding_id_dataset, None)
        if dflt_dim_names != None:
            return dflt_dim_names
        else:
            return SoundingDataFile.get_id_dim_names(self)

    def get_sounding_indexes(self, sounding_id):
        """Find sounding id through bisection, possibly slower than a dict lookup"""
        
        if self.instrument_name == GOSAT_INST_NAME:
            if str(sounding_id)[-1].upper() in GOSAT_POL_ORDER:
                pol_name = str(sounding_id)[-1]
                sounding_id = str(sounding_id)[:-1]
                pol_index = GOSAT_POL_ORDER.index(pol_name)
            else:
                pol_index = slice(len(GOSAT_POL_ORDER))

            index_tuple = SoundingDataFile.get_sounding_indexes(self, sounding_id)
            index_tuple = index_tuple.__class__(index_tuple[0], pol_index)
        else:
            index_tuple = SoundingDataFile.get_sounding_indexes(self, sounding_id)

        return index_tuple

    def determine_instrument_name(self):
        for instrument_name, dataset_name in list(SOUNDING_ID_DATASET.items()):
            if self.get(dataset_name, None) != None:
                return instrument_name
        return None

    def get_info_dataset_name(self, info_name):
        """Looks up the dataset name that matches a given informative short name. The informative name refers
        to a data item that may have a different dataset name between different instruments or might have a
        different name in different versions of the product"""
        
        # Chained to look first in SOUNDING_INFO_DATASETS then ADDL_INFO_DATASET
        item_dataset_spec = SOUNDING_INFO_DATASETS.get(info_name.lower(), ADDL_INFO_DATASETS.get(info_name.lower(), None))

        # Check if the info item undefined or defined for a specific instrument
        if item_dataset_spec == None:
            raise MissingDataset('Could not find dataset for info name: %s' % (info_name))
        elif hasattr(item_dataset_spec, 'get'):
            item_dataset_name = item_dataset_spec.get(self.instrument_name, None)
        else:
            item_dataset_name = item_dataset_spec

        if item_dataset_name == None or len(item_dataset_name) == 0:
            raise MissingDataset('Dataset name for info item: %s and instrument: %s is empty' % (info_name, self.instrument_name))

        return self.select_valid_dataset(item_dataset_name)

    def get_sounding_data(self, data_name, sounding_id=None, indexes=None, flatten=False, average=None, shape_names=None, **kwargs):
        # Get data from base class, snd_data should be a NamedShapeArray class
        snd_data = SoundingDataFile.get_sounding_data(self, data_name, sounding_id, indexes, **kwargs)

        if average != None:
            if average not in snd_data.named_shape._fields:
                raise AveragingError('Can not average over dimension: %s for data named: %s which is not in the dataset indexes: %s' % (average, data_name, snd_data.named_shape))
            else:
                # Average on correct dimension or average the whole array because
                # its shape belongs to the one specified
                if hasattr(snd_data, 'shape'):
                    snd_data = numpy.average(snd_data, snd_data.named_shape._fields.index(average))
                else:
                    snd_data = numpy.average(snd_data)

        if flatten:
            snd_data = numpy.ravel(snd_data)

        return snd_data

                        
    def get_sounding_info(self, info_name, sounding_id=None, **kwargs):

        dataset_name = self.get_info_dataset_name(info_name)
        out_sounding_info = self.get_sounding_data(dataset_name, sounding_id, **kwargs)

        return out_sounding_info

    def get_sounding_info_dict(self, sounding_id, ignore_missing=False, as_strings=False, **kwargs):
        info_dict = {}
        for info_name in list(SOUNDING_INFO_DATASETS.keys()):
            try:

                # If problems averaging, try without it
                try:
                    info_data = self.get_sounding_info(info_name, sounding_id, **kwargs)
                except AveragingError:
                    tmp_kwargs = copy.copy(kwargs)
                    tmp_kwargs['average'] = None
                    info_data = self.get_sounding_info(info_name, sounding_id, **tmp_kwargs)

            except MissingDataset:
                # Ignore missing datasets
                if ignore_missing:
                    continue
                else:
                    raise 

            if as_strings:
                if is_list_like(info_data):
                    info_str = ' '.join([str(value) for value in numpy.ravel(info_data)])
                else:
                    info_str = str(info_data)


                info_dict[info_name] = info_str
            else:
                info_dict[info_name] = info_data


        return info_dict

    def get_sounding_time(self, sounding_id, **kwargs):
        # Calculate current time from TAI93 start
        tmp_kwargs = copy.copy(kwargs)
        tmp_kwargs['average'] = None

        sounding_times = self.get_sounding_data(self.select_valid_dataset(TIME_DATASET), sounding_id, **tmp_kwargs)

        time_structs = []
        for curr_l1b_time in numpy.ravel(sounding_times):
            if isinstance(curr_l1b_time, six.string_types) or type(curr_l1b_time) is numpy.str_:

                if len(curr_l1b_time) < TIME_STRING_EXPECT_LEN:
                    raise Exception('Time string: "%s" from file: "%s" does not have the expected format length: %d' % (curr_l1b_time, self.filename, TIME_STRING_EXPECT_LEN))

                parsed_time = datetime.datetime.strptime(curr_l1b_time[:TIME_STRING_EXPECT_LEN], TIME_STRING_FORMAT)
                unparsed_str = curr_l1b_time[TIME_STRING_EXPECT_LEN:].strip("Z")
                parsed_time = parsed_time + datetime.timedelta(seconds=float(unparsed_str))

                time_structs.append( parsed_time.timetuple() )
            else:
                tai_time = TAI93_START + datetime.timedelta(seconds=(curr_l1b_time + TAI93_LPSECS))
                utc_time = tai64n.tai2utc(tai_time)
                time_structs.append( utc_time.timetuple() )
            
        return tuple(time_structs)

    def get_radiance_data(self, sounding_id, **kwargs):
        radiance_data = []
        for spec_dataset_name in RADIANCE_DATASETS:
            spec_dataset_full = '/%s/%s' % (RADIANCE_GROUP[self.instrument_name], spec_dataset_name)
            
            band_data = self.get_sounding_data(spec_dataset_full, sounding_id, **kwargs)
            
            radiance_data.append( band_data )
                         
        return tuple(radiance_data)

    def get_channel_counts(self, sounding_id, **kwargs):
        return tuple([ len(band_data[0]) for band_data in self.get_radiance_data(sounding_id, **kwargs) ])

    def evaluate_dispersion(self, sounding_id, **kwargs):
        dispersion_coefs = self.get_sounding_info('dispersion', sounding_id, **kwargs)
        channel_counts   = self.get_channel_counts(sounding_id, **kwargs)

        disp_eval = []
        for band_coefs, band_len in zip(dispersion_coefs, channel_counts):
            band_poly = numpy.lib.polynomial.poly1d(band_coefs[::-1])
            disp_eval.append( band_poly(numpy.arange(1,band_len+1)) )

        return tuple(disp_eval)

    def get_wavenumbers(self, sounding_id, **kwargs):
        if self.instrument_name == GOSAT_INST_NAME:
            return self.evaluate_dispersion(sounding_id, **kwargs)
        else:
            return tuple([ old_div(1e4,band_wvl) for band_wvl in self.get_wavelengths(sounding_id, **kwargs) ])

    def get_wavelengths(self, sounding_id, **kwargs):
        if self.instrument_name == OCO_INST_NAME:
            return self.evaluate_dispersion(sounding_id, **kwargs)
        else:
            return tuple([ old_div(1e4,band_wn) for band_wn in self.get_wavenumbers(sounding_id, **kwargs) ])

    def get_error_data(self, sounding_id, calculate_noise=True, gain_code=None, **kwargs):

        # Process noise for each band and return 
        error_data = []
        
        if self.instrument_name == OCO_INST_NAME:
            if gain_code != None:
                raise Exception('gain_code not used for instrument: %s' % OCO_INST_NAME)
            
            index_tuple = self.get_sounding_indexes(sounding_id)
            
            error_data = []
            for band_idx, band_radiance in enumerate(self.get_radiance_data(sounding_id, **kwargs)):
                photon_col = self[OCO_SNR_COEF_DATASET][band_idx, index_tuple[1], :, 0]
                bkgrnd_col = self[OCO_SNR_COEF_DATASET][band_idx, index_tuple[1], :, 1]

                if calculate_noise:
                    tmp = (100.0e0 * band_radiance[:] / MAX_MEAS_SIGNAL[self.instrument_name][band_idx]) * photon_col[:]**2
                    tmp = numpy.sqrt(tmp + bkgrnd_col[:]**2)
                    band_error = (old_div(MAX_MEAS_SIGNAL[self.instrument_name][band_idx],100.0)) * tmp
                else:
                    band_error = (photon_col, bkgrnd_col)

                error_data.append(band_error)
       
        elif self.instrument_name == GOSAT_INST_NAME:

            if gain_code == None:
                gain_code = self.get_sounding_info('gain', sounding_id)
                if is_list_like(gain_code):
                    if not numpy.all(gain_code == gain_code[0]):
                        raise ValueError('sounding id: %s does not specify polarization name and averaging not enababled or gain codes differ for polarization channels: %s' % (sounding_id, gain_code))
                    else:
                        gain_code = gain_code[0]
            
            for band_idx, band_name in enumerate(BAND_DATA_NAMES):
                cnv_dataset = GOSAT_CNV_COEF_DATASET[gain_code.strip()] % band_name
                noise_dataset = GOSAT_NOISE_DATASET % band_name

                cnv_col    = self.get_sounding_data(cnv_dataset, sounding_id, **kwargs)
                band_noise = self.get_sounding_data(noise_dataset, sounding_id, **kwargs)

                if calculate_noise:
                    error_data.append( cnv_col[:] * band_noise )
                else:
                    error_data.append( (cnv_col, band_noise) )
            
        return tuple(error_data)

    def get_build_id(self):
        try:
            b_id_str = self[BUILD_ID_DATASET[self.instrument_name]]
            l1b_build_id = []
            for id_part in b_id_str[0].replace(b'v',b'').split(b'.'):
                if id_part.isdigit():
                    l1b_build_id.append( int(id_part) )
                else:
                    l1b_build_id.append( id_part )
        except KeyError:
            return (0,0,0)
       
        return tuple(l1b_build_id)

    def get_surface_grouping(self, sounding_id):

        surface_type = None

        try:
            # GOSAT and OCO L1B files have differing ways of representing the
            # land type
            if self.instrument_name == GOSAT_INST_NAME:
                land_value = self.get_sounding_info('land_fraction', sounding_id)
                land_values_check_dict = LAND_FRACTION_PERCENTAGE
            elif self.instrument_name == OCO_INST_NAME:
                land_value = self.get_sounding_info('land_water_indicator', sounding_id)
                land_values_check_dict = LAND_WATER_INDICATOR
            else:
                raise MissingDataset("Unknown instrument.")

            for type_name, land_check in list(land_values_check_dict.items()):
                if land_check(land_value):
                    surface_type = type_name
                    break

        except MissingDataset:
            # Fall back to trying to read out of a surface database file
            # Initalize as needed
            if self.surf_type_obj == None:
                self.surf_type_obj = ModisEcoMap()

            # Just use first lat/lon if sounding matches multiple polarizations, etc
            latitude  = self.get_sounding_info('latitude', sounding_id, flatten=True)[0]
            longitude = self.get_sounding_info('longitude', sounding_id, flatten=True)[0]

            surface_type = self.surf_type_obj.get_surface_grouping(latitude, longitude)

        return surface_type


################################################################################################

class L2(SoundingDataFile):
    
    def __init__(self, filename, mode=None, **kwargs):
        # Init super-class
        SoundingDataFile.__init__(self, filename, mode, **kwargs)
        self._sounding_id_dataset = self.select_valid_dataset( ('RetrievalHeader/sounding_id',
                                                                'RetrievalHeader/sounding_id_reference') )

        # These shouldn't be missing from the L2 file
        self._data_shape_name_dict = { "/SpectralParameters/num_colors_per_band": ("Retrieval","Band"),
                                       "/SpectralParameters/measured_radiance": ("Retrieval","Color"),
                                       "/SpectralParameters/measured_radiance_uncert": ("Retrieval","Color"),
                                       "/SpectralParameters/modeled_radiance": ("Retrieval","Color"),
                                       "/SpectralParameters/wavenumber": ("Retrieval","Color"),
                                       "/SpectralParameters/wavelength": ("Retrieval","Color"),
                                       "/RetrievalResults/retrieved_aerosol_aod_by_type": ("Retrieval","AerosolType"),
                                       "/Metadata/AerosolTypes": ("Retrieval", "AerosolName"),
                                       }
        self._default_shape_names = ()

################################################################################################

class IdDatasetFinder(SoundingDataFile):
   
    def __init__(self, filename, mode=None, single_id_dim=False, srch_dim_ids=SOUNDING_ID_DIMENSIONS, **kwargs):
        # Init super-class
        SoundingDataFile.__init__(self, filename, mode, **kwargs)

        self._id_dim_names = None

        if self.get('/SoundingGeometry/sounding_id') != None and len(self['/SoundingGeometry/sounding_id'].shape) == 2:

            # OCO-2 L1B, ECMWF, IMAP  or ABand Cloud screener files
            self._sounding_id_dataset = '/SoundingGeometry/sounding_id'
            if single_id_dim:
                self._id_dim_names = ('Frame',)
            else:
                self._id_dim_names = ('Frame', 'Sounding')

        elif self.get('/RetrievalHeader/sounding_id_reference') != None:
            # L2 output files:
            self._sounding_id_dataset = '/RetrievalHeader/sounding_id_reference'
            self._id_dim_names = ('Retrieval',)

        elif self.get('/RetrievalHeader/sounding_id') != None:
            # L2AggPGE output files:
            self._sounding_id_dataset = '/RetrievalHeader/sounding_id'
            self._id_dim_names = ('Retrieval',)

        elif self.get('/RadianceStatistics_spectra/sounding_id') != None:
            # L1BSt file
            self._sounding_id_dataset = '/RadianceStatistics_spectra/sounding_id'
            if single_id_dim:
                self._id_dim_names = ('Frame',)
            else:
                self._id_dim_names = ('Frame', 'Sounding')

        else:
            for ds_name in list(srch_dim_ids.keys()):
                if self.get(ds_name) != None:
                    self._sounding_id_dataset = ds_name
                    self._id_dim_names = self.get_data_shape(ds_name, srch_dim_ids)

                    if single_id_dim or len(self[ds_name].shape) == 1:
                        self._id_dim_names = (self._id_dim_names[0],)

        if self._id_dim_names == None:
            # Copy filename before closing file. Close file since 
            # this exception means the caller won't have access to
            # this object to close the file
            filename = self.filename
            self.close()
            raise MissingDataset("Unable to find sounding_id dataset for file: %s" % filename)

class SoundingFirstFile(SoundingDataFile):
    def __init__(self, filename, mode=None, single_id_dim=False, **kwargs):
        # Init super-class
        SoundingDataFile.__init__(self, filename, mode, **kwargs)

        # Doesn't exist, instead assume first dim is sounding dimension
        self._sounding_id_dataset = None

        # Specifically return frame only for OCO2 files
        if(self.filename.find("OCO2_") >= 0):
            if single_id_dim:
                self._id_dim_names = ('Frame',)
            else:
                self._id_dim_names = ('Frame', 'Sounding')
        else:
            self._id_dim_names = ('Exposure',) 

    def get_sounding_ids(self):
        # Look for the first dataset in the first group that is not a group itself
        if len(list(self.values())) > 0 and len(list(self.values())[0].values()) > 0 and hasattr(list(list(self.values())[0].values())[0], "shape"):
            indexes = numpy.zeros(list(list(self.values())[0].values())[0].shape, dtype=int)
            indexes.ravel()[:] = list(range(numpy.product(indexes.shape))) 
        else:
            indexes = []
        return indexes 
