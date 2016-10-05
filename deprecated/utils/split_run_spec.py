#!/usr/bin/env python

import os
import sys
import random

from optparse import OptionParser

import L2_Input

if __name__ == "__main__":
    
    # Parse command line options
    parser = OptionParser(usage="usage: %prog [options] <run_spec_file> [ <run_spec_file> ... ]")

    parser.add_option( "-n", "--num_runs", dest="num_runs",
                       metavar='INT',
                       type='int',
                       action='append',
                       help="number of runs per file, can specified more than once"
                       )

    parser.add_option( "-e", "--equal_parts", dest="equal_parts",
                       metavar='INT',
                       type='int',
                       help="split in equal parts per file",
                       )

    parser.add_option( "-o", "--split_names", dest="split_names",
                       metavar='STRING',
                       action='append',
                       help="names to use other than indexes for outputted files"
                       )

    parser.add_option( "-s", "--shuffle", dest="shuffle_runs",
                       default=False,
                       action="store_true",
                       help="Shuffle run directories list before splitting"
                       )

    parser.add_option( "-b", "--binary", dest="binary",
                       help="place supplied binary name into run specs instead of source binary filename",
                       )
    

    # Parse command line arguments
    (cmd_options, cmd_args) = parser.parse_args()

    if cmd_args < 1:
        parser.error('No run_spec file specified')

    if (cmd_options.num_runs == None or len(cmd_options.num_runs) == 0) and cmd_options.equal_parts == None:
        parser.error('Do not know how to split, use -n or -e option')

    if (cmd_options.num_runs != None and len(cmd_options.num_runs) > 0) and cmd_options.equal_parts != None:
        parser.error('Can not specify -n and -e options together')

    input_files = cmd_args

    complete_run_list = []
    headers = []
    for curr_file in input_files:
        if not os.path.exists(curr_file):
            parser.error('run_spec file does not exist: "%s"' % curr_file)

        run_spec_obj = L2_Input.Input_File(curr_file)

        for curr_header in run_spec_obj.Get_Section('HEADER'):
            headers.append(curr_header)

        if run_spec_obj.Get_Matrix_Data() != None:
            [ complete_run_list.append(run_dir) for run_dir in run_spec_obj.Get_Matrix_Data() ]

    if cmd_options.shuffle_runs:
        random.shuffle(complete_run_list)            

    if cmd_options.equal_parts != None:
        num_runs = [ cmd_options.equal_parts for index in range(len(complete_run_list) / cmd_options.equal_parts + 1) ]
    else:    
        num_runs = cmd_options.num_runs
        
    num_runs_idx = 0
    file_run_lists = []
    while len(complete_run_list) > 0:
        if num_runs_idx < len(num_runs):
            num_to_pull = num_runs[num_runs_idx]
            if len(complete_run_list) > num_to_pull:
                
                file_run_lists.append( complete_run_list[0:num_to_pull] )
                complete_run_list = complete_run_list[num_to_pull:]

            else:
                file_run_lists.append( complete_run_list )    
                complete_run_list = []
            num_runs_idx += 1
        else:
            file_run_lists.append( complete_run_list )    
            complete_run_list = []
            
    print 'Creating %d new run_spec files' % len(file_run_lists)

    idx_padding = len(str(len(file_run_lists)))

    out_uniq_name = os.path.commonprefix([ os.path.basename(in_name) for in_name in input_files ])

    dot_pos = out_uniq_name.rfind('.')
    if dot_pos < len(out_uniq_name) - 5:
        new_ext = '.dat'
        dot_pos = len(out_uniq_name) - len(new_ext)
        out_uniq_name += new_ext
    
    for run_list_idx in range(len(file_run_lists)):
        if cmd_options.split_names != None and run_list_idx < len(cmd_options.split_names):
            postfix = cmd_options.split_names[run_list_idx]
        else:
            postfix = ('%0' + str(idx_padding) + 'd') % (run_list_idx)

        output_filename = out_uniq_name[0:dot_pos] + '_' + postfix + out_uniq_name[dot_pos:]
        file_runs = file_run_lists[run_list_idx]
        
        print 'Creating "%s" with %d runs' % (output_filename, len(file_runs))

        new_run_spec = L2_Input.Input_File()

        new_header = L2_Input.Section('HEADER')

        new_header.Set_Keyword_Value('num_runs', str(len(file_runs)))
        new_header.Set_Keyword_Value('source_run_spec', ', '.join([ os.path.realpath(in_file) for in_file in input_files ]))
            
        for head_sect in headers:
            for key_name in head_sect.Get_All_Keyword_Names():
                new_header.Set_Keyword_Value(key_name, head_sect.Get_Keyword_Value(key_name))

        if cmd_options.binary != None:
            new_header.Set_Keyword_Value('binary_filename', os.path.realpath(cmd_options.binary))

        new_run_spec.children.append(new_header)

        new_run_spec.Set_Matrix_Data(file_runs)
        new_run_spec.Write(output_filename, doIndent=True)
