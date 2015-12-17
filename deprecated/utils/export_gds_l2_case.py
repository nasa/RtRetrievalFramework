#!/usr/bin/env python

import sys
import os
from optparse import OptionParser

import L2_Input

def get_external_paths(file_obj):
    sections = file_obj.rootNode.Get_All_Sections()

    print sections

path_input_files = ['oco_l2.inp', 'oco_l2.run']

parser = OptionParser(usage="usage: %prog [options] <testcase_directory> [testcase_directory]...")

# Parse command line arguments
(cmd_options, testcase_dirs) = parser.parse_args()

if len(testcase_dirs) == 0:
    parser.error('at least one testcase directory must be specified')


for tc_dir in testcase_dirs:
    tc_dir = os.path.realpath(tc_dir)

    for inp_file_base in path_input_files:

        inp_file_name = tc_dir + "/" + inp_file_base

        if not os.path.exists(inp_file_name):
            print '%s does not exist for testcase dir %s' % (inp_file_base, tc_dir)
            continue

        inp_file_obj = L2_Input.Input_File(inp_file_name)

        get_external_paths(inp_file_obj)

