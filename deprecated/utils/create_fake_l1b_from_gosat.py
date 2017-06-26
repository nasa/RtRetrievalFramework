#!/usr/bin/env python

import os
import sys
import time
import contextlib
from optparse import OptionParser

import numpy
import h5py

import GOSAT_File

CONVERT_GUIDE = { 'FootprintGeometry': [ { 'dataset': 'footprint_latitude',
                                           'type'   : numpy.float32,
                                           'query'  : GOSAT_File.L1B.get_swir_geometry,
                                           'keyword': 'latitude'
                                           },
                                         { 'dataset': 'footprint_longitude',
                                           'type'   : numpy.float32,
                                           'query'  : GOSAT_File.L1B.get_swir_geometry,
                                           'keyword': 'longitude'
                                           },
                                         { 'dataset': 'footprint_solar_zenith',
                                           'type'   : numpy.float64,
                                           'query'  : GOSAT_File.L1B.get_swir_geometry,
                                           'keyword': 'solar_zenith',
                                           'add_dim': 2,
                                           },

                                         { 'dataset': 'footprint_time',
                                           'type'   : numpy.float64,
                                           'query'  : GOSAT_File.L1B.get_oco_times,
                                           },

                                         { 'dataset': 'footprint_altitude',
                                           'type'   : numpy.float32,
                                           }
                                         ],
                  'FrameHeader'      : [ { 'dataset': 'frame_id',
                                           'type'   : numpy.int64,
                                           'query'  : GOSAT_File.L1B.get_observation_id,
                                           },
                                         { 'dataset': 'frame_time_stamp',
                                           'type'   : '|S25',
                                           'query'  : GOSAT_File.L1B.get_time_stamp,
                                           'add_dim': 1,
                                           },
                                         ],
                  'SoundingGeometry' : [ { 'dataset': 'sounding_id',
                                           'type'   : numpy.int64,
                                           'query'  : GOSAT_File.L1B.get_observation_id,
                                           'add_dim': 1,
                                           }
                                         ],
                  }

def create_fake_l1b(gosat_l1b_file, oco_l1b_file):

    with contextlib.nested( contextlib.closing(GOSAT_File.L1B(gosat_l1b_file)),
                            contextlib.closing(h5py.File(oco_l1b_file, 'w')) ) as (gosat_l1b_obj, oco_l1b_obj):

        # Set up reported soundings
        metadata = oco_l1b_obj.require_group('Metadata')
        reported_snds = numpy.ones(1, dtype=numpy.int8)
        metadata.create_dataset('ReportedSoundings', data=reported_snds)

        num_obs = gosat_l1b_obj.get_num_observations()
        for group_name, convert_items in CONVERT_GUIDE.items():
            base_group = oco_l1b_obj.require_group(group_name)
            print 'Created group: %s' % base_group
        
            for convert_settings in convert_items:
                if convert_settings.has_key('query'):
                    
                    for obs_index in range(num_obs):
                        print 'Copying %s at %s' % (convert_settings['dataset'],
                                                    gosat_l1b_obj.get_time_stamp(obs_index))
                        
                        in_data = convert_settings['query'](gosat_l1b_obj, obs_index)

                        if convert_settings.has_key('keyword'):
                            in_data = in_data[ convert_settings['keyword'] ]

                        if obs_index == 0:                           
                            new_dims = [num_obs] + list( getattr(in_data, 'shape', []) )
                            if convert_settings.has_key('add_dim'):
                                for dim_idx in range(convert_settings['add_dim']):
                                    new_dims.append(1)
                            new_data = numpy.zeros(new_dims, dtype=convert_settings['type'])

                        if len(new_dims)-1 == 0:
                            new_data[obs_index] = in_data
                        elif len(new_dims)-1 == 1:
                            new_data[obs_index, :] = in_data
                        elif len(new_dims)-1 == 2:
                            new_data[obs_index, :, :] = in_data
                        else:
                            raise Exception('%s not created due to unsupported shape of input data: %s' % (convert_settings['dataset'], str(new_dims)))

                new_dataset = base_group.create_dataset(convert_settings['dataset'], data=new_data)
                print 'Created dataset: %s' % str(new_dataset)


        
def standalone_main():
    parser = OptionParser(usage="usage: %prog [options] hdf_filename [output_file]")

    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('Need to at least specify output file')

    gosat_l1b_file = args[0]

    if len(args) < 2:
        oco_l1b_file = os.path.splitext(gosat_l1b_file)[0] + '_OCO.hdf'
        print 'Writing output to: %s' % oco_l1b_file
    else:
        oco_l1b_file = args[1]

    create_fake_l1b(gosat_l1b_file, oco_l1b_file)

if __name__ == "__main__":
    standalone_main()
