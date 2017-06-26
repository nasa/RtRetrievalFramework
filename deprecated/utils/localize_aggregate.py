#!/usr/bin/env python

import os
import re
import sys
import shutil

from optparse import OptionParser

from L2_FP_Util import run_wait_command
import L2_Input

NEW_CONFIG_FILENAME    = '%s.local'
DEFAULT_REMOTE_MACHINE = 'graphite.jpl.nasa.gov'
INPUT_FILE_SECTIONS    = [ 'input->InputProductFiles',
                           'input->InputProductFiles->SpectrumFiles'
                           ]
USED_FILENAME_KEYWORDS = [ 'L1BFile',
                           'ResampledMetFile',
                           'RunlogFile',
                           '\d{17}'
                           ]

SSH_EXEC = 'ssh -q -o ForwardX11=no'

if __name__ == "__main__":
    
    # Parse command line options
    parser = OptionParser(usage="usage: %prog [options] <agg_config_file> [ <agg_config_file> ... ]")

    parser.add_option( "-e", "--external_input_dir", dest="external_input_dir",
                       metavar='DIR',
                       default='./external_input',
                       help='directory other than current where to copy external input files',
                       )

    parser.add_option( "-a", "--absolute_paths", dest="absolute_paths",
                       default=False,
                       action="store_true",
                       help="replace path to external input files in config files with absolute rather than relative paths"
                       )

    parser.add_option( "-r", "--remote_machine", dest="remote_machine",
                       metavar='HOSTNAME',
                       default=DEFAULT_REMOTE_MACHINE,
                       help='alternative remote machine where aggregate files live',
                       )


    # Parse command line arguments
    (cmd_options, cmd_args) = parser.parse_args()

    if len(cmd_args) == 0:
        parser.error('No aggregate config filenames specified')

    external_inp_dir = os.path.realpath(cmd_options.external_input_dir)

    # Make sure all the files exist before processing any
    for agg_config_file in cmd_args:
        if not os.path.exists(agg_config_file):
            parser.error('The aggregate config file: %s does not exist' % agg_config_file)

    # List of files to scp at one time
    external_copy_list = []

    # Save aggregate objects only after all error checks passed
    agg_objs_to_save = []
    
    for agg_config_file in cmd_args:
        agg_dirname = os.path.dirname(os.path.realpath(agg_config_file))

        # Open as an object to read filenames and change programmatically
        print 'Processing config file: %s' % agg_config_file
        agg_obj = L2_Input.Input_File(agg_config_file)

        agg_objs_to_save.append(agg_obj)

        for inp_sect_name in INPUT_FILE_SECTIONS:
            for inp_sect_obj in agg_obj.Get_Section(inp_sect_name):
                for inp_key_name in inp_sect_obj.Get_All_Keyword_Names():

                    # Ignore files we do not need to copy
                    is_key_used = False
                    for key_check in USED_FILENAME_KEYWORDS:
                        if re.search(key_check, inp_key_name):
                            is_key_used = True
                            break
                        
                    if not is_key_used:
                        continue
                    
                    ext_fullname = inp_sect_obj.Get_Keyword_Value(inp_key_name)

                    if ext_fullname == None or len(ext_fullname) == 0:
                        raise IOError('Keyword %s in section: %s of file: "%s" is empty' % (inp_key_name, inp_sect_name, agg_config_file))

                    # Only change things in the config file if the file does not exist locally
                    if not os.path.exists(ext_fullname):
                        if ext_fullname not in external_copy_list:
                            external_copy_list.append(ext_fullname)

                        file_basename = os.path.basename(ext_fullname)

                        if cmd_options.absolute_paths:
                            local_fullname = os.path.join(external_inp_dir, file_basename)
                        else:
                            local_fullname = os.path.join(os.path.relpath(external_inp_dir, agg_dirname), file_basename)

                        print '%s:' % inp_key_name
                        print '%s -> %s' % (ext_fullname, local_fullname)
                        print ''
                        inp_sect_obj.Set_Keyword_Value(inp_key_name, local_fullname)

    if len(external_copy_list) == 0:
        print >>sys.stderr, 'No valid files found in config files for copying'
        sys.exit(1)
        
    # Check that the files to copy all exist on the foreign system
    test_expressions = ' -a '.join([ '"(" -e "%s" ")"' % remote_file for remote_file in external_copy_list ])
    test_command = 'ssh -q %s test \'%s\'' % (cmd_options.remote_machine, test_expressions)

    print 'Testing for file existance on machine: %s' % cmd_options.remote_machine
    error_code = run_wait_command(test_command, False)
    
    if error_code != 0:
        print >>sys.stderr, 'One or more of the following remote files do not exist on machine: %s' % cmd_options.remote_machine
        for remote_file in external_copy_list:
            print >>sys.stderr, remote_file
        print >>sys.stderr, 'One or more of the previous remote files do not exist on machine: %s' % cmd_options.remote_machine
        sys.exit(error_code)
    else:
        print 'Verified existance of %d remote files' % len(external_copy_list)

    # Make destination directory if it does not exist
    if not os.path.exists(external_inp_dir):
        os.makedirs(external_inp_dir)

    # Copy files from remote system
    copy_list_str = ' '.join([ '"%s:%s"' % (cmd_options.remote_machine, remote_file) for remote_file in external_copy_list ])
    copy_command = 'rsync -Cauzbv --rsh "%s" %s %s' % (SSH_EXEC, copy_list_str, external_inp_dir)

    print 'Copying %d files to: %s' % (len(external_copy_list), external_inp_dir)
    error_code = run_wait_command(copy_command, verbose=False, show_output=True)
    
    if error_code != 0:
        print >> sys.stderr, 'Error occurred copying %d files from remote machine: %s' % (len(external_copy_list), cmd_options.remote_machine)

    # Verify that all files were copied
    all_copied = True
    for remote_filename in external_copy_list:
        file_basename = os.path.basename(remote_filename)
        local_filename = os.path.join(external_inp_dir, file_basename)
        if not os.path.exists(local_filename):
            print >>sys.stderr, '%s was not succesfully copied' % local_filename
            all_copied = False

    # If no errors have occured yet so update config files
    if all_copied:
        print 'Saving updated aggregate config files'
        for agg_obj in agg_objs_to_save:
            agg_obj.Write(NEW_CONFIG_FILENAME % agg_obj.filename, doIndent=True)
    
