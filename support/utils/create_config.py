#!/usr/bin/env python

import os
import re
import sys
import traceback

from optparse import OptionParser

import h5py

from full_physics import L2InputFile, parse_keyval_str_list, PopulatorBase
import full_physics.acos_file as acos_file


FRAME_ID_DATASET = { acos_file.OCO_INST_NAME: '/FrameHeader/frame_id',
                     acos_file.GOSAT_INST_NAME: '/SoundingHeader/sounding_id',
                     }
SOUNDING_ID_DATASET = { acos_file.OCO_INST_NAME: '/SoundingGeometry/sounding_id',
                        acos_file.GOSAT_INST_NAME: '/SoundingHeader/sounding_id',
                        }

LATITUDE_DATASET  = 'FootprintGeometry/footprint_latitude'
LONGITUDE_DATASET = 'FootprintGeometry/footprint_longitude'

SZA_DATASET = 'FootprintGeometry/footprint_solar_zenith'

# Group names for detection that do not appear in the other products
L1B_SRCH_GROUP = 'InstrumentHeader'
ECMWF_SRCH_GROUPS = ['ecmwf', 'ECMWF']
IMAP_SRCH_GROUP = 'DOASFluorescence'
CLOUD_SRCH_GROUP = 'ABandCloudScreen'
RRV_SRCH_GROUP = 'campaign'

# What to modify in input config templates
INP_PROD_SECTION_NAME     = 'input->InputProductFiles'
PGE_GROUP_SECTION_NAME    = 'input->PGENameGroup'
ALT_SETTINGS_SECTION_NAME = 'input->AlternativeSettings'

GEN_LIST_SECTION_TMPL     = 'input->%sFullPhysics'
SOUNDING_ID_LIST_NAME     = 'SoundingIds'

def insert_alternative_settings(template_obj, alt_settings_hash):

    alt_sect = template_obj.get_section(ALT_SETTINGS_SECTION_NAME)

    if len(alt_sect) == 0:
        raise IOError('Could not section %s inside of template file %s' % (ALT_SETTINGS_SECTION_NAME, template_obj.filename))
    else:
        for (alt_key, alt_val) in alt_settings_hash.items():
            alt_sect[0].set_keyword_value(alt_key, alt_val)

def find_valid_file(check_list, ext_check_dict):

    for curr_file in check_list:
        (file_prefix, file_ext) = os.path.splitext(curr_file)
        if ext_check_dict.has_key(file_ext[1:]): # 1: removes . at beg of extension
            check_result = ext_check_dict[file_ext[1:]](curr_file)
            if check_result != None:
                return check_result 

    return None

def process_l1b_sounding_ids(template_obj, l1b_obj, filter_funcs=[], l1b_keyword='L1BFile', id_list_sect=None):
    
    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)

    if len(inp_prod_section) == 0:
        print template_obj.get_all_section_names()
        raise IOError('Could not find input product file section of %s' % template_obj.filename)

    inp_prod_section[0].set_keyword_value(l1b_keyword, l1b_obj.filename)

    snd_id_matrix = l1b_obj.get_sounding_ids()
    num_snd_dims = len(snd_id_matrix.shape)

    config_ids = []
    for row_idx in range(snd_id_matrix.shape[0]):
        if num_snd_dims == 2:
            num_cols = snd_id_matrix.shape[1]
        else:
            num_cols = 1
        for col_idx in range(num_cols):
            if num_snd_dims == 2:
                curr_sounding_ids  = [ snd_id_matrix[row_idx, col_idx] ]
                curr_sounding_idxs = [ (row_idx, col_idx) ]
            else:
                curr_sounding_ids  = [ snd_id_matrix[row_idx] ]
                curr_sounding_idxs = [ (row_idx,) ]

            for curr_filter in filter_funcs:
                curr_sounding_ids,curr_sounding_idxs = curr_filter(curr_sounding_ids, curr_sounding_idxs)
       
            if curr_sounding_ids != None:
                config_ids += curr_sounding_ids

    sounding_ids_val_sect = []
    if id_list_sect == None:
        sounding_id_list_sect = GEN_LIST_SECTION_TMPL % l1b_obj.instrument_name
    else:
        sounding_id_list_sect = id_list_sect

    for list_sect in template_obj.get_section(sounding_id_list_sect + '->LIST'):
        list_name = list_sect.get_keyword_value('name')
        if list_name == SOUNDING_ID_LIST_NAME:
            sounding_ids_val_sect = list_sect.get_section('->VALUES')

    if len(sounding_ids_val_sect) == 0:
        raise IOError('Could not find sounding id list section named %s in %s' % (sounding_id_list_sect, template_obj.filename))

    sounding_ids_val_sect[0].set_matrix_data(config_ids)


def handle_common_config(template_obj, out_config_filename, used_files, filter_options, sounding_filters=[]):

    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)

    if len(inp_prod_section) == 0:
        print template_obj.get_all_section_names()
        raise IOError('Could not find input product file section of %s' % template_obj.filename)

    def check_l1b_file(check_file):
        with h5py.File(check_file, 'r') as hdf_obj:
            if L1B_SRCH_GROUP in  hdf_obj.keys():
                return check_file
            else:
                return None

    def check_ecmwf_file(check_file):
        with h5py.File(check_file, 'r') as hdf_obj:
            for ecmwf_grp in ECMWF_SRCH_GROUPS:
                if ecmwf_grp  in hdf_obj.keys():
                    return check_file
            return None

    def check_imap_file(check_file):
        with h5py.File(check_file, 'r') as hdf_obj:
            if IMAP_SRCH_GROUP in hdf_obj.keys():
                return check_file
            else:
                return None

    def check_cloud_file(check_file):
        with h5py.File(check_file, 'r') as hdf_obj:
            if CLOUD_SRCH_GROUP in hdf_obj.keys():
                return check_file
            else:
                return None

    def check_rrv_file(check_file):
        with h5py.File(check_file, 'r') as hdf_obj:
            if RRV_SRCH_GROUP in hdf_obj.keys():
                return check_file
            else:
                return None

    # L1B files come either with .hdf or .h5 extension depending
    # on when it was processed
    l1b_file = find_valid_file(used_files, { 'hdf': check_l1b_file, 'h5': check_l1b_file })
    ecmwf_file = find_valid_file(used_files, { 'hdf': check_ecmwf_file, 'h5': check_ecmwf_file })
    imap_file = find_valid_file(used_files, { 'hdf': check_imap_file, 'h5': check_imap_file })
    cloud_file = find_valid_file(used_files, { 'hdf': check_cloud_file, 'h5': check_cloud_file })
    rrv_file = find_valid_file(used_files, { 'hdf': check_rrv_file, 'h5': check_rrv_file })

    if l1b_file != None:
        l1b_obj = acos_file.L1B(l1b_file)

        def sounding_arg_filters(in_sounding_ids, in_sounding_idxs):

            out_sounding_ids = []
            for curr_sounding_id in in_sounding_ids:
                valid_sounding = True
                curr_pos = str(curr_sounding_id)[-1]

                if filter_options['sounding_positions'] != None and len(filter_options['sounding_positions']) != 0 and not curr_pos in filter_options['sounding_positions']:
                    valid_sounding = False

                if filter_options['sounding_ids'] != None and len(filter_options['sounding_ids']) != 0 and not str(curr_sounding_id) in filter_options['sounding_ids']:
                    valid_sounding = False

                if valid_sounding:
                    out_sounding_ids.append(curr_sounding_id)
                    
            return out_sounding_ids, in_sounding_idxs

        def surface_grouping_filter(in_sounding_ids, in_sounding_idxs):
            out_sounding_ids  = []
            out_sounding_idxs = []
            for curr_sounding_id, curr_sounding_idxs in zip(in_sounding_ids, in_sounding_idxs):

                try:
                    surf_grouping = l1b_obj.get_surface_grouping(curr_sounding_id)
                except ValueError:
                    print >>sys.stderr, 'Error retrieving surface type for %s' % curr_sounding_id
                    traceback.print_exception(*sys.exc_info(), limit=2)
                    break

                for curr_surf_group in filter_options['surface_groupings']:
                    if curr_surf_group == surf_grouping:
                        print curr_sounding_id, 'recognized as surface grouping:', surf_grouping
                        out_sounding_ids.append(curr_sounding_id)
                        out_sounding_idxs.append(curr_sounding_idxs)
                        break

            return out_sounding_ids, out_sounding_idxs

        def solar_zenith_angle_filter(in_sounding_ids, in_sounding_idxs):
            out_sounding_ids  = []
            out_sounding_idxs = []
            for curr_sounding_id, curr_sounding_idxs in zip(in_sounding_ids, in_sounding_idxs):
                if len(curr_sounding_idxs) == 3:
                    sza = (l1b_obj[SZA_DATASET][curr_sounding_idxs[0], curr_sounding_idxs[1], curr_sounding_idxs[2]],)
                else:
                    sza  = l1b_obj[SZA_DATASET][curr_sounding_idxs[0], :, :].ravel()

                if sza[0] <= int(filter_options['solar_zenith_angle']):
                    print curr_sounding_id, 'has sza %d <= %d' % (sza[0], int(filter_options['solar_zenith_angle']))
                    out_sounding_ids.append(curr_sounding_id)
                    out_sounding_idxs.append(curr_sounding_idxs)

            return out_sounding_ids, in_sounding_idxs


        def cloud_flag_filter(in_sounding_ids, in_sounding_idxs):
            out_sounding_ids  = []
            out_sounding_idxs = []

            cloud_obj = h5py.File(cloud_file, "r")
            cloud_flag = cloud_obj['/ABandCloudScreen/cloud_flag']
            for curr_sounding_id, curr_sounding_idxs in zip(in_sounding_ids, in_sounding_idxs):
                if cloud_flag.__getitem__(curr_sounding_idxs) == 0:
                    out_sounding_ids.append(curr_sounding_id)
                    out_sounding_idxs.append(curr_sounding_idxs)

            return out_sounding_ids, in_sounding_idxs

        if filter_options.get('sounding_positions') != None or filter_options.get('sounding_ids') != None:
            sounding_filters.append(sounding_arg_filters)

        if filter_options.get('surface_groupings') != None:
            sounding_filters.append(surface_grouping_filter)

        if filter_options.get('solar_zenith_angle') != None:
            sounding_filters.append(solar_zenith_angle_filter)

        if filter_options.get('cloud_flag') and cloud_file != None:
            sounding_filters.append(cloud_flag_filter)

        process_l1b_sounding_ids(template_obj, l1b_obj, sounding_filters)

        # Only close after processing because filter functions use the file
        l1b_obj.close()
           
    if ecmwf_file != None:
        inp_prod_section[0].set_keyword_value('ResampledMetFile', ecmwf_file)

    if imap_file != None:
        inp_prod_section[0].set_keyword_value('IMAPFile', imap_file)

    if cloud_file != None:
        inp_prod_section[0].set_keyword_value('CloudFile', cloud_file)

    if rrv_file != None:
        inp_prod_section[0].set_keyword_value('RRVFile', rrv_file)

    template_obj.write(out_config_filename, doIndent=True)

# Same as gen_config with defaults
handle_oco_config = handle_common_config

def handle_oco_uplooking_config(template_obj, out_config_filename, used_files, filter_options):

    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)[0]                                                                                                        
    for curr_file in used_files:
        extension = curr_file[curr_file.rfind('.')+1:].lower()
        
        if extension == 'dat':
            inp_prod_section.set_keyword_value('AtmosphereFile', curr_file)

    sounding_filters = []
    handle_common_config(template_obj, out_config_filename, used_files, filter_options, sounding_filters)

def handle_gosat_config(template_obj, out_config_filename, used_files, filter_options):

    def filter_pol_posistion(in_sounding_ids, in_sounding_idxs):
        out_sounding_ids = []
        out_sounding_idxs = []
        for curr_sounding_id, curr_sounding_idx in zip(in_sounding_ids, in_sounding_idxs):
            id_string = '%d' % curr_sounding_id
            out_sounding_ids.append( id_string + 'P' )
            out_sounding_ids.append( id_string + 'S' )

            out_sounding_idxs.append( (curr_sounding_idx[0], 0, 0) )
            out_sounding_idxs.append( (curr_sounding_idx[0], 0, 1) )
            
        return out_sounding_ids, out_sounding_idxs

    sounding_filters = []
    if (filter_options['sounding_ids'] == None or len(filter_options['sounding_ids']) == 0) and filter_options['separate_pol']:
        sounding_filters.append(filter_pol_posistion)
    handle_common_config(template_obj, out_config_filename, used_files, filter_options, sounding_filters)

def handle_fts_config(template_obj, out_config_filename, used_files, filter_options):

    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)[0] 
    spectrum_files_section = inp_prod_section.get_section('->SpectrumFiles')[0]
    
    obs_ids = []
    for curr_file in used_files:
        extension = curr_file[curr_file.rfind('.')+1:].lower()
        
        if extension == 'grl':
            inp_prod_section.set_keyword_value('RunlogFile', curr_file)
        if extension == 'dat':
            inp_prod_section.set_keyword_value('AtmosphereFile', curr_file)
        if re.search('\d\d\d', extension):
            if filter_options['sounding_ids'] != None and not extension in filter_options['sounding_ids']:
                continue
            
            if not extension in obs_ids:
                obs_ids.append(extension)
            
            curr_id = extension
            if re.search('^.*a[_.]', curr_file):
                curr_id += '1'
            elif re.search('^.*b[_.]', curr_file):
                curr_id += '2'

            spectrum_files_section.set_keyword_value(curr_id, curr_file)

    obs_ids.sort()
    obs_id_sec = template_obj.get_section('input->FTSFullPhysics->LIST->VALUES')[0]
    obs_id_sec.set_matrix_data(obs_ids)

    template_obj.write(out_config_filename, doIndent=True)

def handle_uq_config(template_obj, out_config_filename, used_files, filter_options):

    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)[0] 

    uq_obj = acos_file.L1B(used_files[0])

    process_l1b_sounding_ids(template_obj, uq_obj, l1b_keyword='UqFile', id_list_sect='input->UqFullPhysics')

    template_obj.write(out_config_filename, doIndent=True)

if __name__ == "__main__":
    
    # Parse command line options
    parser = OptionParser(usage="usage: %prog [options] -t <run_type> <type_file> [<type_file>...]")

    parser.add_option( "-t", "--run_type", dest="run_type",
                       metavar='STRING', 
                       help="type of run being set up, valid values are: [%s]" % (', '.join(PopulatorBase.populator_list.keys())),
                       )

    parser.add_option( "-l", "--sounding_list", dest="sounding_list_file",
                       metavar='FILE',
                       help="sounding list to use when picking sounding ids placed into config file"
                       )


    parser.add_option( "-s", "--sounding_pos", dest="sounding_positions",
                       metavar='STRING',
                       action='append',
                       help="filter sounding ids by the specified sounding posistion, specify multiple times for multiple sounding files"
                       )

    parser.add_option( "-g", "--surface_groupings", dest="surface_groupings",
                       metavar='STRING',
                       action='append',
                       help="filter sounding ids by the surface grouping returned from Modis ECO Map"
                       )

    parser.add_option( "-c", "--cloud_flag", dest="cloud_flag",
                       default=False,
                       action="store_true",
                       help="Use cloud file (if detected) cloud flag to filter soundings")

    parser.add_option( "-z", "--sza_cutoff", dest="sza_cutoff",
                       metavar='INT',
                       help="filter sounding ids for solar zenith angles less than cutoff"
                       )

    parser.add_option( "-p", "--separate_pol", dest="separate_pol",
                       default=False,
                       action="store_true",
                       help="for GOSAT sounding ids, split ids from L1B into P and S ids ")

    parser.add_option( "-o", "--config_filename", dest="config_filename",
                       metavar='FILE',
                       help="alternative name for output config file produced other than based on run type and directory name"
                       )

    parser.add_option( "-a", "--alternative_setting", dest="alternative_settings",
                       metavar="KEY=VALUE",
                       type="string",
                       action="append",
                       help='define an override control setting keyword and value. Multiple are allowed and are defined using the syntax: "keyword=value"'
                       )

    parser.add_option( "-r", "--relative_paths", dest="relative_paths",
                       default=False,
                       action="store_true",
                       help='Store paths in config file as relative to directory where config file will end up'
                       )

    # Parse command line arguments
    (cmd_options, cmd_args) = parser.parse_args()

    if len(cmd_args) < 1:
        parser.error('No files specified for setting up config')

    if cmd_options.run_type == None:
        parser.error('run type not set with -t argument')

    if not PopulatorBase.populator_list.has_key(cmd_options.run_type):
        parser.error('%s is not a valid run type, valid values are: [%s]' % (cmd_options.run_type, ', '.join(PopulatorBase.populator_list.keys())))

    if cmd_options.config_filename == None:
        out_config_filename = '%s/%s_%s.config' % (os.getcwd(), cmd_options.run_type, os.path.basename(os.getcwd()))
    else:
        out_config_filename = cmd_options.config_filename

    if cmd_options.relative_paths:
        cfg_path = os.path.dirname(out_config_filename)
        used_files = [ os.path.relpath(curr_file, cfg_path) for curr_file in cmd_args ]
    else:
        used_files = [ os.path.realpath(curr_file) for curr_file in cmd_args ]

    populator = PopulatorBase.create_populator_from_config_type(cmd_options.run_type)
    if not os.path.exists(populator.config_template_file):
        raise IOError('Template config file does not exist: %s' % populator.config_template_file)
    
    tmpl_obj = L2InputFile(populator.config_template_file)

    sounding_ids = None
    if cmd_options.sounding_list_file != None:
        sounding_ids = populator.read_id_list_file(cmd_options.sounding_list_file)

    filter_options = { 'sounding_positions': cmd_options.sounding_positions,
                       'sounding_ids':       sounding_ids,
                       'surface_groupings':  cmd_options.surface_groupings,
                       'solar_zenith_angle': cmd_options.sza_cutoff,
                       'separate_pol':       cmd_options.separate_pol,
                       'cloud_flag':         cmd_options.cloud_flag,
                       }

    # Insert any alternative configuration values
    alt_settings_hash = parse_keyval_str_list(cmd_options.alternative_settings)
    insert_alternative_settings(tmpl_obj, alt_settings_hash)

    eval('handle_%s_config(tmpl_obj, out_config_filename, used_files, filter_options)' % cmd_options.run_type)
