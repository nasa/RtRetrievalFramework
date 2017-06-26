from __future__ import division
from builtins import range
from builtins import object
from past.utils import old_div

# Reads CERES/SARB surface properties maps from:
# http://www-surf.larc.nasa.gov/surf/pages/data-page.html
#
# "All maps, except the digital elevation map are 8-bit binary data made
# on an SGI workstation. Their size is 2160 points in longitude, 1080
# points in latitude or 1/6 degree equal angle. All maps begin at the
# North Pole, Greenwich Meridian."

import os
import struct
import math

import array
import bisect
import numpy
import h5py


RECORD_FORMAT = '>c'

LATITUDE_RANGE  = [-90, 90]
LONGITUDE_RANGE = [-180, 180]

FRACTION_37 = 0.5
PERIOD = math.pi*2.0/365.0
    
LONGITUDE_CENTER = 0

HDF_MAP_TYPE = 'modis'

class SurfaceBinMap(object):
    """Base class for handling common read/write operations from
    CERES/SARB surface property maps"""

    map_resolution = { 'modis': 60.0,
                       'ceres':  6.0,
                       }
    
    def __init__(self, filename):

        self.filename = filename

        if not os.path.exists(self.filename):
            raise IOError('%s map file does not exist' % self.filename)
        
        # Determine map type and parameters
        self._determine_map_params()
        
        # Open file for access
        self._open()


    def __del__(self):

        self.close()

    def _open(self):

        self.map_f_obj = open(self.filename, 'rb')       

    def close(self):

        self.map_f_obj.close()

    def _determine_map_params(self):

        # Get file info, size is in bytes
        stat_info = os.stat(self.filename)

        self.map_type = None
        for map_type, resolution in list(self.map_resolution.items()):
            self._NUM_LONGITUDE = (LONGITUDE_RANGE[1] - LONGITUDE_RANGE[0]) * resolution
            self._NUM_LATITUDE  = (LATITUDE_RANGE[1]  - LATITUDE_RANGE[0])  * resolution
            self._DELTA_ANGLE   = resolution
            curr_type_size = self._NUM_LONGITUDE * self._NUM_LATITUDE

            if curr_type_size == stat_info.st_size:
                self.map_type = map_type
                break
         
        if self.map_type == None:
            raise IOError('Could not determine map type for file: %s' % (self.filename))
        

    def _check_lat_lon(self, latitude, longitude):
        if not (latitude >= LATITUDE_RANGE[0] and latitude <= LATITUDE_RANGE[1]):
            raise ValueError('Latitude: %f out of range: [%f, %f]' % (latitude, LATITUDE_RANGE[0], LATITUDE_RANGE[1]))
           
        if not (longitude >= LONGITUDE_RANGE[0] and longitude <= LONGITUDE_RANGE[1]):
            raise ValueError('Longitude: %f out of range: [%f, %f]' % (longitude, LONGITUDE_RANGE[0], LONGITUDE_RANGE[1]))

    def _read_value(self, lon_idx, lat_idx):

        if self.map_type == 'ceres':
            record_idx = int( lat_idx*self._NUM_LONGITUDE+lon_idx )
        elif self.map_type == 'modis':
            record_idx = int( lon_idx*self._NUM_LATITUDE+lat_idx )
            #record_idx = int( lat_idx*self._NUM_LONGITUDE+lon_idx )

        self.map_f_obj.seek(record_idx,0)

        # Read binary data into an integer from file
        rec_str  = self.map_f_obj.read(struct.calcsize(RECORD_FORMAT))
        rec_data = struct.unpack(RECORD_FORMAT, rec_str)
        data_value = ord(rec_data[0])

        return data_value

    def _check_indexes(self, lat_idx, lon_idx, latitude, longitude):

        if not (lat_idx >= 0 and lat_idx < self._NUM_LATITUDE):
            raise ValueError('%d is an invalid latitude index for latitude %f' % (lat_idx, latitude))

        if not (lon_idx >= 0 and lon_idx < self._NUM_LONGITUDE):
            raise ValueError('%d is an invalid longitude index for longitude %f' % (lon_idx, longitude))

    def get_value(self, latitude, longitude, surrounding=None, **kwds):

        # Quick check of lat/lon values within accepted ranges
        self._check_lat_lon(latitude, longitude)

        # Convert latitude/longitude to x,y
        if self.map_type == 'ceres':
            # Map left is greenich
            if longitude < LONGITUDE_CENTER:
                i_lon = int( ((LONGITUDE_RANGE[1] - LONGITUDE_RANGE[0]) + longitude) * self._DELTA_ANGLE )
            else:
                i_lon = int( longitude * self._DELTA_ANGLE )

            # Map top is north pole
            i_lat = int( ((LATITUDE_RANGE[1] - LATITUDE_RANGE[0]) - (latitude  - LATITUDE_RANGE[0])) * self._DELTA_ANGLE )
        elif self.map_type == 'modis':
            # Map left is -180
            i_lon = round( (longitude - (old_div(2.0,self._DELTA_ANGLE) + LONGITUDE_RANGE[0])) * self._DELTA_ANGLE )
            i_lat = round( (latitude  - (old_div(2.0,self._DELTA_ANGLE) + LATITUDE_RANGE[0])) * self._DELTA_ANGLE )
            
        else:
            raise ValueError('%s is an unknown map type' % self.map_type)

        self._check_indexes(i_lat, i_lon, latitude, longitude)

        if surrounding != None:
            delta_range = list(range(-surrounding,surrounding+1))

            values = numpy.zeros((len(delta_range),len(delta_range)), dtype=int)
            for lon_delta in delta_range:
                for lat_delta in delta_range:
                    values[lat_delta-delta_range[0], lon_delta-delta_range[0]] = self._read_value(i_lon+lon_delta, i_lat+lat_delta, *kwds)
            return values
        else:
            return self._read_value(i_lon, i_lat, *kwds)

class IgbpTypeMap(object):
    """Handles IGBP+1 CERES scene type maps"""

    type_to_name = {
         1: 'Evergreen Needleleaf Forest',
         2: 'Evergreen Broadleaf Forest',
         3: 'Deciduous Needleleaf Forest',
         4: 'Deciduous Broadleaf Forest',
         5: 'Mixed Forest', # 'Mixed Deciduous Forest'
         6: 'Closed Shrublands',
         7: 'Open Shrubland(Desert)',
         8: 'Woody Savanna',
         9: 'Savanna',
        10: 'Grassland',
        11: 'Permenant Wetland',
        12: 'Cropland',
        13: 'Urban',
        14: 'Crop Mosaic', # 'Crop/Natural Veg: "Mosaic"'
        15: 'Permanent Snow',
        16: 'Barren/Desert',
        17: 'Water',
        18: 'Tundra',
        19: 'Fresh Snow',
        20: 'Sea Ice',
        }

    def __init__(self, filename):
        # Create map of the types and names in reverse fashion
        self.name_to_type = {}
        for type_num, name_str in list(self.type_to_name.items()):
            self.name_to_type[name_str] = type_num


class ModisEcoMap(IgbpTypeMap):

    MAX_SURF_MASK_VALUE = 7

    mask_to_name = {
        0: 'Shallow ocean',
        1: 'Land',
        2: 'Ocean coastlines and lake shorelines',
        3: 'Shallow inland water',
        4: 'Ephemeral water',
        5: 'Deep inland water',
        6: 'Moderate or continental ocean',
        7: 'Deep ocean',
        }  

    surface_groupings = {
        'land':  [ 'Land',
                   'Ocean coastlines and lake shorelines',                  
                   'Shallow inland water',
                   'Ephemeral water',
                   'Deep inland water',
                   ],
        'ocean': [ 'Shallow ocean',
                   'Moderate or continental ocean',
                   'Deep ocean', ]
        }

    def __init__(self, filename=None):
        if filename == None:
            if not 'L2_DATASETS_PATH' in os.environ:
                raise Exception("L2_DATASETS_PATH not set, please see user's guide for setup instructions")

            filename = os.path.join(os.environ['L2_DATASETS_PATH'], 'surface/modis/IGBP.EcoMap.v1.0.2004.129.v004.h5')

        if not os.path.exists(filename):
            raise IOError('Can not load non existant map file: %s' % filename)
        
        # Add type mappings specific to this class
        self.type_to_name[0]   = 'Water'
        self.type_to_name[254] = 'Unclassified'

        # Create reverse of mask_to_name to map name to a mask
        self.name_to_mask = {}
        for type_num, name_str in list(self.mask_to_name.items()):
            self.name_to_mask[name_str] = type_num

        # Create mapping of mask to grouping name
        self.mask_to_grouping = {}
        for grouping_name, grouping_types in list(self.surface_groupings.items()):
            for mask_id, mask_name in list(self.mask_to_name.items()):
                if mask_name in grouping_types:
                    self.mask_to_grouping[mask_id] = grouping_name

        # Check that we mapped all types to groupings
        if len(self.mask_to_grouping) != len(self.mask_to_name):
            raise Exception('Could not map all mask types to surface grouping names')

        IgbpTypeMap.__init__(self, filename)

        self.filename = filename
        self.map_f_obj = h5py.File(self.filename, 'r')

        self.latitudes_rev = list(self.map_f_obj['Latitude'])
        self.latitudes_rev.reverse()

    def close(self):
        self.map_f_obj.close()

    def get_value(self, latitude, longitude, dataset='IGBP_Land_Cover_Type'):

        # Latitude in HDF file are sorted in decreasing order so
        # latitudes_rev is reversed from the HDF file into sorted order
        # Then we use bisect to find the index in this list and then
        # reverse the index back to one that matches the HDF file.
        # Got that? I had forgotten it too before I added this comment.
        lat_idx_fnd = bisect.bisect_left(self.latitudes_rev, latitude)
        lat_idx_fnd = len(self.latitudes_rev) - lat_idx_fnd

        # Thankfully longitudes are sorted in increasing order
        lon_idx_fnd = bisect.bisect_left(self.map_f_obj['Longitude'], longitude)

        # Bisect does not always give us the closest lat/lon indexes so look around
        # 2 points on either side and see if there are any closer matches.
        lat_idx_min = max(0, min(lat_idx_fnd, self.map_f_obj['Latitude'].shape[0]-1))
        lon_idx_min = max(0, min(lon_idx_fnd, self.map_f_obj['Longitude'].shape[0]-1))
        for offset in (-1, 1):
            
            lat_idx_curr = lat_idx_fnd+offset
            if 0 <= lat_idx_curr < self.map_f_obj['Latitude'].shape[0]:
                lat_diff_curr = abs(latitude  - self.map_f_obj['Latitude'][lat_idx_curr])
                lat_diff_min = abs(latitude  - self.map_f_obj['Latitude'][lat_idx_min])
            
                if lat_diff_curr < lat_diff_min:
                    lat_idx_min = lat_idx_curr

            lon_idx_curr = lon_idx_fnd+offset
            if 0 <= lon_idx_curr < self.map_f_obj['Longitude'].shape[0]:
                lon_diff_curr = abs(longitude - self.map_f_obj['Longitude'][lon_idx_curr])
                lon_diff_min = abs(longitude - self.map_f_obj['Longitude'][lon_idx_min])

                if lon_diff_curr < lon_diff_min:
                    lon_idx_min = lon_idx_curr

        return self._read_value(lon_idx_min, lat_idx_min, dataset)
       
    def _read_value(self, lon_idx, lat_idx, dataset):

        return self.map_f_obj[dataset][lat_idx][lon_idx]

    def get_igbp_type(self, latitude, longitude):

        return self.get_value(latitude, longitude, dataset='IGBP_Land_Cover_Type')

    def get_land_cover_dict(self, latitude, longitude):

        cover_value = self.get_value(latitude, longitude, dataset='Land_Cover_Type_QC')

        cover_dict = {}
        cover_dict['byte_value']      = cover_value
        cover_dict['mandatory_qa']    = cover_value >> 0 & 3  # first two bits
        cover_dict['quarters_since']  = cover_value >> 2 & 3  # next two bits
        cover_dict['land_water_mask'] = cover_value >> 4 & 15  # last four bits

        if cover_dict['land_water_mask'] > self.MAX_SURF_MASK_VALUE:
            raise ValueError('Extracted wrong land water mask value: %d for latitude: %s, longitude: %s from file: %s' % (cover_dict['land_water_mask'], latitude, longitude, self.filename))
            
        return cover_dict

    def get_surface_grouping(self, latitude, longitude):

        land_water_mask = self.get_land_cover_dict(latitude, longitude)['land_water_mask']
        
        return self.mask_to_grouping[land_water_mask]
        
#######
    
    # Routines to retrieve IGBP scene classification, to return albedo over
    # land/water (get_albedo) and to return albedo over snow/ice
    # (get_snow_ice_albedo)
    #
    # Routines as implemented by Eric Moody and Steve Platnick.  See original
    # code for details

class AlbedoTable(object):
    
    table_wavelengths = (0.659, 0.858, 1.24, 1.64, 2.13, 3.74)

    #     Load the albedo tables
    #     Define the summer albedo values for the ecosystem types:
    #     0.659     0.858    1.24     1.64     2.13     3.74 microns
    albedo_summer = numpy.array([
        [ 0.04222, 0.22031, 0.20521, 0.10548, 0.04508, 0.04508 ], #Evergreen Needle Forest
        [ 0.05603, 0.33227, 0.32429, 0.19531, 0.08681, 0.08681 ], #Evergreen Broad  Forest
        [ 0.08036, 0.23529, 0.25979, 0.20509, 0.12427, 0.12427 ], #Deciduous Needle Forest
        [ 0.10243, 0.26507, 0.28033, 0.22043, 0.12955, 0.12955 ], #Deciduous Broad  Forest
        [ 0.07133, 0.27526, 0.26623, 0.15418, 0.07286, 0.07286 ], #Mixed Forest
        [ 0.08412, 0.24195, 0.27441, 0.23508, 0.14127, 0.14127 ], #Closed Shrubland
        [ 0.14669, 0.24542, 0.30695, 0.31364, 0.24182, 0.24182 ], #Open Shrubland
        [ 0.07045, 0.28070, 0.30356, 0.22453, 0.12145, 0.12145 ], #Woody Savannas
        [ 0.08328, 0.28865, 0.32313, 0.25405, 0.14490, 0.14490 ], #Savannas
        [ 0.10992, 0.27952, 0.33016, 0.29120, 0.18843, 0.18843 ], #Grasslands
        [ 0.07998, 0.20398, 0.20124, 0.14967, 0.08590, 0.08590 ], #Perm. wetlands
        [ 0.11641, 0.30108, 0.34447, 0.29252, 0.17672, 0.17672 ], #Croplands
        [ 0.09087, 0.25731, 0.26656, 0.20708, 0.13707, 0.13707 ], #Urban & built-up
        [ 0.07417, 0.34395, 0.35212, 0.24735, 0.12635, 0.12635 ], #Cropland mosaics
        [ 0.90000, 0.85000, 0.45000, 0.15000, 0.15000, 0.05000 ], #Snow/Ice
        [ 0.34885, 0.41770, 0.48565, 0.52709, 0.50776, 0.25388 ], #Barren/Sparse Veg.
        [ 0.05000, 0.05000, 0.05000, 0.05000, 0.05000, 0.05000 ], #Water Bodies
        [ 0.10000, 0.25000, 0.30000, 0.25000, 0.15000, 0.05000 ], #Tundra
        ])
    
    #     Define the winter albedo values for the ecosystem types:
    #     0.659     0.858    1.24     1.64     2.13     3.74 microns
    albedo_winter = numpy.array([
        [ 0.04549, 0.17819, 0.17644, 0.10573, 0.05247, 0.05247 ], #Evergreen Needle Forest
        [ 0.04552, 0.31785, 0.31662, 0.18155, 0.07232, 0.07232 ], #Evergreen Broad  Forest
        [ 0.09192, 0.18807, 0.19381, 0.13663, 0.08652, 0.08652 ], #Deciduous Needle Forest
        [ 0.05761, 0.23056, 0.26090, 0.18531, 0.09542, 0.09542 ], #Deciduous Broad  Forest
        [ 0.05678, 0.19229, 0.21550, 0.15289, 0.08270, 0.08270 ], #Mixed Forest
        [ 0.08755, 0.21747, 0.25394, 0.20686, 0.13120, 0.13120 ], #Closed Shrubland
        [ 0.16818, 0.24631, 0.29580, 0.29723, 0.24372, 0.24372 ], #Open Shrubland
        [ 0.07508, 0.23402, 0.26787, 0.20593, 0.12080, 0.12080 ], #Woody Savannas
        [ 0.09578, 0.23720, 0.27900, 0.23363, 0.14931, 0.14931 ], #Savannas
        [ 0.13212, 0.21712, 0.28075, 0.28429, 0.21291, 0.21291 ], #Grasslands
        [ 0.07509, 0.21504, 0.22970, 0.15472, 0.08098, 0.08098 ], #Perm. wetlands
        [ 0.09031, 0.22343, 0.26643, 0.23283, 0.15232, 0.15232 ], #Croplands
        [ 0.09776, 0.21394, 0.23763, 0.20206, 0.13869, 0.13869 ], #Urban & built-up
        [ 0.07863, 0.25835, 0.29249, 0.22286, 0.12651, 0.12651 ], #Cropland mosaics
        [ 0.90000, 0.85000, 0.45000, 0.15000, 0.15000, 0.05000 ], #Snow/Ice
        [ 0.34885, 0.41770, 0.48565, 0.52709, 0.50776, 0.25388 ], #Barren/Sparse Veg.
        [ 0.05000, 0.05000, 0.05000, 0.05000, 0.05000, 0.05000 ], #Water Bodies
        [ 0.10000, 0.25000, 0.30000, 0.25000, 0.15000, 0.05000 ], #Tundra
        ])
    
    #     Define the tropical albedo values for the ecosystem types:
    #     0.659     0.858    1.24     1.64     2.13     3.74 microns
    albedo_tropical = numpy.array([
        [ 0.05999, 0.22086, 0.22479, 0.14854, 0.07236, 0.07236 ], #Evergreen Needle Forest
        [ 0.06772, 0.34919, 0.34498, 0.19382, 0.08220, 0.08220 ], #Evergreen Broad  Forest
        [ 0.08141, 0.22551, 0.25907, 0.20851, 0.12899, 0.12899 ], #Deciduous Needle Forest
        [ 0.10812, 0.28096, 0.34392, 0.28955, 0.18194, 0.18194 ], #Deciduous Broad  Forest
        [ 0.05836, 0.25317, 0.27197, 0.18123, 0.08823, 0.08823 ], #Mixed Forest
        [ 0.10439, 0.27004, 0.32049, 0.27573, 0.17893, 0.17893 ], #Closed Shrubland
        [ 0.18282, 0.29743, 0.36292, 0.35973, 0.28930, 0.28930 ], #Open Shrubland
        [ 0.07166, 0.28618, 0.31932, 0.22378, 0.11542, 0.11542 ], #Woody Savannas
        [ 0.09004, 0.27721, 0.32717, 0.26052, 0.15267, 0.15267 ], #Savannas
        [ 0.12031, 0.26763, 0.33087, 0.29940, 0.20617, 0.20617 ], #Grasslands
        [ 0.09627, 0.29647, 0.29990, 0.19487, 0.10267, 0.10267 ], #Perm. wetlands
        [ 0.08767, 0.31173, 0.33267, 0.24248, 0.13566, 0.13566 ], #Croplands
        [ 0.09662, 0.25265, 0.26485, 0.20728, 0.13839, 0.13839 ], #Urban & built-up
        [ 0.08569, 0.31492, 0.34767, 0.25772, 0.14099, 0.14099 ], #Cropland mosaics
        [ 0.90000, 0.85000, 0.45000, 0.15000, 0.15000, 0.05000 ], #Snow/Ice
        [ 0.34885, 0.41770, 0.48565, 0.52709, 0.50776, 0.25388 ], #Barren/Sparse Veg.
        [ 0.05000, 0.05000, 0.05000, 0.05000, 0.05000, 0.05000 ], #Water Bodies
        [ 0.10000, 0.25000, 0.30000, 0.25000, 0.15000, 0.05000 ], #Tundra
        ])
    
    # Set the 3.7 value to be a fraction of the 2.1 micron values:
    albedo_summer[6,1:14]   = albedo_summer[5,1:14] * FRACTION_37
    albedo_winter[6,1:14]   = albedo_winter[5,1:14] * FRACTION_37
    albedo_tropical[6,1:14] = albedo_tropical[5,1:14] * FRACTION_37
    
    albedo_dry_snow = numpy.array( (0.900, 0.850, 0.450, 0.150, 0.150, 0.050) )
    albedo_wet_snow = numpy.array( (0.850, 0.750, 0.250, 0.050, 0.050, 0.050) )

    def get_albedo(self, jday, latitude, longitude, wavelengths, scene_type):

        scene_type = int(scene_type)

        jday = min(365, jday)
              
        if scene_type < 19:
            #  Test to see if this is in the tropics, defined as between +- 20
            #  degrees.
            if ( latitude >= -20.0 and latitude <= 20.0 ):
                #  This pixel is in the tropics.  In the tropics, we do not have
                #  a seasonal variation, so just set the albedo to the
                #  corresponding ecosystem type.
                albedos = self.albedo_tropical[:, scene_type]

            elif ( latitude > 20.0 and latitude < 30.0 ):
                # This pixel is in a transition zone between northern hemisphere
                # and tropical. To avoid a sharp discontinuity, a linear fit is
                # used between the northern hemisphere and tropical albedo values.
                # First the northern hemisphere's albedo is computed, and then
                # the linear fit is performed.  For details of the northern
                # hemisphere albedo computation, see the last else section.

                # Define the summer solstice date:
                summer_solstice = 173

                # Compute the Amplitude, y_offset, hemisphere albedo,
                # slope and intercept, and store the computed albedo value.

                # Define the Amplitude:
                amplitude = old_div((self.albedo_summer[scene_type,:] -
                             self.albedo_winter[scene_type,:]), 2.0)

                # Define the y_offset:
                y_offset = self.albedo_summer[scene_type,:] - amplitude[:]

                # Store the computed albedo value:
                alb_hem = amplitude[:] * \
                          math.cos(PERIOD * (float(jday) - float(summer_solstice))) + \
                            y_offset[:]

                # Compute the slope:
                slope = old_div((alb_hem[:] - self.albedo_tropical[scene_type,:]), 10.0)

                # Compute the intercept:
                intercept = 3.0 * self.albedo_tropical[scene_type,:] - 2.0 * alb_hem[:]

                # Compute and store the albedo:
                albedos = slope[:] * latitude + intercept[:]

            elif ( latitude < -20.0 and latitude > -30.0 ):
                # This pixel is in a transition zone between southern hemisphere and tropical.  To
                # avoid a sharp discontinuity, a linear fit is used between southern
                # hemisphere and tropical albedo values. First the southern hemisphere's albedo
                # is computed, and then the linear fit is performed.  For details of the southern
                # hemisphere albedo computation, see the last else section.

                # Define the summer solstice date:
                summer_solstice = 356

                # Compute the Amplitude, y_offset, hemisphere albedo,
                # slope and intercept, and store the computed albedo value.

                # Define the Amplitude:
                amplitude = old_div((self.albedo_summer[scene_type,:] - \
                             self.albedo_winter[scene_type,:]), 2.0)

                # Define the y_offset:
                y_offset = self.albedo_summber[scene_type,:] - amplitude[:]

                # Store the computed albedo value:
                alb_hem = amplitude[:] * \
                          math.cos(PERIOD*(float(jday)-float(summer_solstice))) + \
                          y_offset[:]

                # Compute the slope:
                slope = old_div((self.albedo_tropical[scene_type,:] - alb_hem[:]), 10.0)

                # Compute the intercept:
                intercept = 3.0*self.albedo_tropical[scene_type,:] - 2.0*alb_hem[:]

                # Compute and store the albedo:
                albedos = slope[:] * latitude + intercept[:]

            else:
                # This pixel resides either in the northern or southern hemisphere,
                # defined as either > 30 or < -30 degrees latitude.  There are
                # representative values for winter and summer albedo values, and
                # to create a seasonal variability, a cosine fit will be used.
                # The fit is based upon the julian day for the summer solstice
                # (June 22st, day 173, for the Northern Hemisphere, and December
                # 22nd, day 356, for the Southern Hemisphere). The form of the fit
                # is as follows:
                # Albedo = Amplitude * cos(Period * (jday - Summer_solstice))
                #          + y_offset
                # where:
                # Period = 2*pi / number of days in year,
                # number of days in year = 365
                # Summer solstice = 172 (NH), or 356 (SH)
                # Amplitude = (summer solstice alb val -
                #                        winter solstice alb val) / 2
                # y_offset = summer solstice alb val - Amplitude

                # Define the Summer solstice date:
                if (latitude >= 30.0):
                    # Northern Hemisphere
                    summer_solstice = 173
                else:
                    # Southern Hemisphere
                    summer_solstice = 356
                    
                # Compute the Amplitude, y_offset, and store the
                # computed albedo value.

                # Define the Amplitude:
                amplitude = old_div((self.albedo_summer[scene_type,:] - \
                             self.albedo_winter[scene_type,:]), 2.0)
                # Define the y_offset:
                y_offset = self.albedo_summer[scene_type,:] - \
                           amplitude[:]

                # Store the computed albedo value:
                albedos = amplitude[:] * \
                          math.cos(PERIOD*(float(jday)-float(summer_solstice))) + \
                          y_offset[:]
        else:
            #  Snow or ice
            albedos = self.compute_si_albedo(jday,latitude,scene_type)

        return self.interp_wavelengths(albedos, wavelengths)

    def compute_si_albedo(self, doy, lat, snow_ice_type):

        if 152 <= doy <= 244:
            arctic_melting_season = True
        else:
            arctic_melting_season = False

        if lat < 0.0: arctic_melting_season = not arctic_melting_season

        albedo_out = None
        if snow_ice_type == 19: #dry snow
            albedo_out = self.albedo_dry_snow[:]
        elif snow_ice_type == 20:   #sea ice
            if arctic_melting_season:
                albedo_out = self.albedo_wet_snow[:]
            else:
                albedo_out = self.albedo_dry_snow[:]
        elif snow_ice_type == 21:   #permanent ice
            albedo_out = self.albedo_dry_snow[:]
        elif snow_ice_type == 22:   #wet snow
            albedo_out = self.albedo_wet_snow[:]

        return albedo_out

    # Perform a linear interpolation in wavelength to get the
    # output albedos at the desired wavelengths.
    def interp_wavelengths(self, albedos, wavelengths):

        # Get albedo at each output wavelength.
        albedo_out = []
        for curr_wl in wavelengths:
            where_wl = numpy.where(curr_wl <= self.table_wavelengths)
            if len(where_wl) == 0:
                wl_table_idx = len(self.table_wavelengths)-2
            else:
                wl_table_idx = min(self.table_wavelengths[0][0], len(self.table_wavelengths)-2)

            denom = self.table_wavelengths[wl_table_idx+1]-self.table_wavelengths[wl_table_idx]

            if abs(denom) < 1e-5:
                frac = 0.0
            else:
                frac = old_div((curr_wl - self.table_wavelengths[wl_table_idx]), denom)

            albedo_out.append( (1.0-frac)*albedos[wl_table_idx]+frac*albedos[wl_table_idx+1] )

        return numpy.array( albedo_out )
                               

