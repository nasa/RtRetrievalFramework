#!/usr/bin/env python

from __future__ import print_function
import os
import re

from optparse import OptionParser

import h5py

from full_physics import L2InputFile, parse_keyval_str_list, PopulatorBase
import full_physics.acos_file as acos_file

# Maps the output locations in the config with the groups in the HDF files indentifying the file type
FILE_SEARCH_GROUPS = { "L1BFile" :   ['InstrumentHeader'],
                       "ResampledMetFile" : ['ecmwf', 'ECMWF',
                                             'Meteorology/meteorology_flag'],
                       "CO2PriorFile" : ['CO2Prior'],
                       "IMAPFile" :  ['DOASFluorescence'],
                       "ABandFile" : ['ABandRetrieval'],
                       "RRVFile" :   ['campaign'],
                       'UqFile'  :   ['StateVector'],
                      }

HDF_VALID_EXTENSIONS = ['.h5', '.hdf']

INPUT_FILE_MAP_EXT = ".map"
INPUT_FILE_MAP_KEYWORD = "InputFileMapping"

# What to modify in input config templates
INP_PROD_SECTION_NAME     = 'input->InputProductFiles'
ALT_SETTINGS_SECTION_NAME = 'input->AlternativeSettings'

GEN_LIST_SECTION_TMPL     = 'input->%sFullPhysics'
SOUNDING_ID_LIST_NAME     = 'SoundingIds'

def insert_alternative_settings(template_obj, alt_settings_hash):

    alt_sect = template_obj.get_section(ALT_SETTINGS_SECTION_NAME)

    if len(alt_sect) == 0:
        raise IOError('Could not section %s inside of template file %s' % (ALT_SETTINGS_SECTION_NAME, template_obj.filename))
    else:
        for (alt_key, alt_val) in list(alt_settings_hash.items()):
            alt_sect[0].set_keyword_value(alt_key, alt_val)

def check_file_type(filename):
    (file_prefix, file_ext) = os.path.splitext(filename)

    file_type_keywords = []
    if file_ext in HDF_VALID_EXTENSIONS:
        # A file could match multiple types, for instance if several types of information are wrapped in one file
        with h5py.File(filename, 'r') as hdf_obj:
            for keyword_name, groups in list(FILE_SEARCH_GROUPS.items()):
                for curr_group in groups:
                    if curr_group in hdf_obj:
                        file_type_keywords.append(keyword_name)

    elif file_ext == INPUT_FILE_MAP_EXT:
        file_type_keywords.append(INPUT_FILE_MAP_KEYWORD)
    else:
        return None

    return file_type_keywords

def handle_common_config(template_obj, out_config_filename, used_files, sounding_ids, run_type, ids_file_keyword='L1BFile', id_list_sect=None):
    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)

    if len(inp_prod_section) == 0:
        print(template_obj.get_all_section_names())
        raise IOError('Could not find input product file section of %s' % template_obj.filename)

    for curr_file in used_files:
        keyword_names = check_file_type(curr_file)
        if keyword_names is not None:
            for kw in keyword_names:
                inp_prod_section[0].set_keyword_value(kw, curr_file)
    
    ids_file = inp_prod_section[0].get_keyword_value(ids_file_keyword)
    if ids_file != None and len(ids_file) > 0 and ids_file != "NONE" and sounding_ids == None:
        l1b_obj = acos_file.L1B(ids_file)
        sounding_ids = list(l1b_obj.get_sounding_ids()[:].ravel())

    sounding_ids_val_sect = []

    if id_list_sect == None:
        id_list_sect = GEN_LIST_SECTION_TMPL % run_type.upper()

    for list_sect in template_obj.get_section(id_list_sect + '->LIST'):
        list_name = list_sect.get_keyword_value('name')
        if list_name == SOUNDING_ID_LIST_NAME:
            sounding_ids_val_sect = list_sect.get_section('->VALUES')

    if len(sounding_ids_val_sect) == 0:
        raise IOError('Could not find sounding id list section named %s in %s' % (id_list_sect, template_obj.filename))

    sounding_ids_val_sect[0].set_matrix_data(sounding_ids)
   
    template_obj.write(out_config_filename, doIndent=True)

# Same as gen_config with defaults
handle_oco_config = handle_common_config
handle_oco3_config = handle_common_config

def handle_oco_uplooking_config(template_obj, out_config_filename, used_files, sounding_ids, run_type):

    inp_prod_section = template_obj.get_section(INP_PROD_SECTION_NAME)[0]                                                                                                        
    for curr_file in used_files:
        extension = curr_file[curr_file.rfind('.')+1:].lower()
        
        if extension == 'dat':
            inp_prod_section.set_keyword_value('AtmosphereFile', curr_file)

    handle_common_config(template_obj, out_config_filename, used_files, sounding_ids, run_type, id_list_sect='input->OCOFullPhysics')

def handle_gosat_config(template_obj, out_config_filename, used_files, sounding_ids, run_type):

    handle_common_config(template_obj, out_config_filename, used_files, sounding_ids, run_type)

def handle_fts_config(template_obj, out_config_filename, used_files, sounding_ids, run_type):

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
            if sounding_ids != None and not extension in sounding_ids:
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

def handle_uq_config(template_obj, out_config_filename, used_files, filter_options, run_type):

    handle_common_config(template_obj, out_config_filename, used_files, sounding_ids, run_type, ids_file_keyword='UqFile', id_list_sect='input->UqFullPhysics')

if __name__ == "__main__":
    
    # Parse command line options
    parser = OptionParser(usage="usage: %prog [options] -t <run_type> <type_file> [<type_file>...]")

    parser.add_option( "-t", "--run_type", dest="run_type",
                       metavar='STRING', 
                       help="type of run being set up, valid values are: [%s]" % (', '.join(list(PopulatorBase.populator_list.keys()))),
                       )

    parser.add_option( "-l", "--sounding_list", dest="sounding_list_file",
                       metavar='FILE',
                       help="sounding list to use when picking sounding ids placed into config file"
                       )

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

    if cmd_options.run_type not in PopulatorBase.populator_list:
        parser.error('%s is not a valid run type, valid values are: [%s]' % (cmd_options.run_type, ', '.join(list(PopulatorBase.populator_list.keys()))))

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

    # Insert any alternative configuration values
    alt_settings_hash = parse_keyval_str_list(cmd_options.alternative_settings)
    insert_alternative_settings(tmpl_obj, alt_settings_hash)

    eval('handle_%s_config(tmpl_obj, out_config_filename, used_files, sounding_ids, cmd_options.run_type)' % cmd_options.run_type)
