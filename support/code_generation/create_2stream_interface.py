#!/usr/bin/env python

import os

import fparser.api

from fp_interface_tools import *

##

import argparse

parser = argparse.ArgumentParser(description='Create 2stream C++ interfaces')
parser.add_argument('twostream_path', help='Path to 2steam source code')
parser.add_argument('output_path', help='Path where output files will be written')
args = parser.parse_args()

##

TWOSTREAM_MAIN_DIR = args.twostream_path
MASTER_FILE_MODULES = {
    os.path.join(TWOSTREAM_MAIN_DIR, '2stream_lps_master.f90'): 'twostream_lps_master_m',
    os.path.join(TWOSTREAM_MAIN_DIR, '2stream_ls_brdf_supplement.f90'): 'twostream_ls_brdf_supplement_m',
    }

MASTERS_IGNORE_ROUTINES = ( 'twostream_lps_fourier_master', 
                            'twostream_ls_brdf_maker',
                            'twostream_ls_brdf_function',
                            'twostream_ls_brdf_fourier',
                            )

SIZE_VARIABLE_NAMES = [ "nlayers",
                        "ntotal",
                        "nbeams",
                        "n_user_streams",
                        "n_user_relazms",
                        "n_geometries",
                        "nstreams_brdf",
                        "maxtotal",
                        "maxmessages",
                        "max_brdf_kernels",
                        "max_brdf_parameters",
                        "maxstreams_brdf",
                        "maxbeams",
                        "max_user_streams",
                        "max_user_obsgeoms",
                        "max_atmoswfs",
                        "max_surfacewfs",
                        "max_user_relazms",
                        "maxlayers",
                        "max_sleavewfs",
                        "max_geometries",
                        ]

MASTERS_CONSTRUCTOR_ARGUMENTS={}
MASTERS_CONSTRUCTOR_ARGUMENTS["twostream_lps_master_m"] = [ 'earth_radius' ] + SIZE_VARIABLE_NAMES
MASTERS_CONSTRUCTOR_ARGUMENTS["twostream_ls_brdf_supplement_m"] = SIZE_VARIABLE_NAMES

NON_ATTRIBUTE_TYPES={}

ACCESSOR_CONST = {}
ACCESSOR_CONST["twostream_l_master_m"] = lambda x: not x.find("brdf_f") >= 0

##
# Output control

INTERFACE_MASTERS_NAME = "twostream_interface"
F_INTERFACE_MASTERS_FILENAME = os.path.join(args.output_path, "%s.F90" % INTERFACE_MASTERS_NAME)
H_INTERFACE_MASTERS_FILENAME = os.path.join(args.output_path, "%s.h" % INTERFACE_MASTERS_NAME)
I_INTERFACE_MASTERS_FILENAME = os.path.join(args.output_path, "%s.i" % INTERFACE_MASTERS_NAME)
TST_INTERFACE_MASTERS_FILENAME = os.path.join(args.output_path, "%s_test.cc" % INTERFACE_MASTERS_NAME)

#### 

# Parse master files
def twostream_variable_configuration(routine_obj, var_obj):
    var_conf = { # Use name itself in fortran code since the 
        # size for arrays are passed along with the array
        "f_arg_var_name": "{name}",
        }

    if var_obj.name in SIZE_VARIABLE_NAMES:
        var_conf["size_variable"] = True
    else:
        var_conf["size_variable"] = False
    
    if var_obj.name in MASTERS_CONSTRUCTOR_ARGUMENTS.get(routine_obj.parent.name, []):
        var_conf["constructor_arg"] = True
        var_conf["class_attribute"] = True
    elif type(var_obj.typedecl) is ftypes.Type and not var_obj.name in NON_ATTRIBUTE_TYPES.get(routine_obj.parent.name, []):
        var_conf["class_attribute"] = True
    elif not var_obj.name in NON_ATTRIBUTE_TYPES.get(routine_obj.parent.name, []):
        var_conf["class_attribute"] = True

    # Base shape and len variables name itself instead of arg name
    if var_conf["class_attribute"]:
        var_conf["c_shape_var_name"] = "{name}_shape_{{0}}"
        var_conf["c_len_var_name"] = "{name}_len"

    # Make non-const accessors for certain variables
    is_accessor_const = ACCESSOR_CONST.get(routine_obj.parent.name, lambda x: True)
    var_conf["accessor_return_const"] = is_accessor_const(var_obj.name)
    var_conf["accessor_method_const"] = is_accessor_const(var_obj.name)

    # Use a different variable for the lower bound of these string arrays
    if var_obj.name.lower() == "c_messages" or var_obj.name.lower() == "c_actions":
        var_conf["f_shape_loop_var_name"] = "min(c_nmessages, {f_shape_var_name})"

    return var_conf

masters_list = parse_master_files(MASTER_FILE_MODULES,
                                  variable_configuration=twostream_variable_configuration,
                                  ignore_routines=MASTERS_IGNORE_ROUTINES)

###
# Create C masters classes
write_cpp_master_classes(H_INTERFACE_MASTERS_FILENAME, I_INTERFACE_MASTERS_FILENAME, masters_list, has_read_write=False)

###
# Create F masters wrapper
write_fortran_master_modules(F_INTERFACE_MASTERS_FILENAME, masters_list, addl_modules=[])
