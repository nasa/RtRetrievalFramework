#!/usr/bin/env python

import os
import sys

# Requires F2PY G3 code from http://code.google.com/p/f2py/
# Used changeset:   81:69b8613d30cb
from fparser import api as fapi, block_statements as block

from fp_interface_tools import *

try:
    # Python 2.7+
    from collections import OrderedDict
except ImportError:
    # http://pypi.python.org/pypi/ordereddict
    from ordereddict import OrderedDict

##
# Fortran input filenames

LIDORT_BASE_DIR = os.path.join(os.environ["L2_EXE_PATH"], "thirdparty/lidort-3.6/")
LIDORT_DEFS_DIR = os.path.join(LIDORT_BASE_DIR, "lidort_def")
LIDORT_MAIN_DIR = os.path.join(LIDORT_BASE_DIR, "lidort_main")
BRDF_MAIN_DIR = os.path.join(LIDORT_BASE_DIR, "sup/brdf")

PARS_FILENAME = "lidort_pars.F90.in"
LIDORT_PARS_NAME = "lidort_pars"
LIDORT_INTERFACE_TYPES_IO_NAME = "lidort_interface_types_io"

DEF_BASE_CLASS_NAME = "Lidort_Structure"

# Filenames with files containing type definitions 
DEF_FILENAMES = []

# BRDF type definitions
BRDF_DEF_FILENAMES = (
    "brdf_lin_sup_inputs_def.F90",
    "brdf_lin_sup_outputs_def.F90",
    "brdf_sup_inputs_def.F90",
    "brdf_sup_outputs_def.F90",
        )

for brdf_def_fn in BRDF_DEF_FILENAMES:
    DEF_FILENAMES.append(os.path.join(BRDF_MAIN_DIR, brdf_def_fn))

# Derived type definitions
LIDORT_DEF_FILENAMES = (
    "lidort_io_defs.F90",
    "lidort_lin_inputs_def.F90",
    "lidort_lin_io_defs.F90",
    "lidort_lin_outputs_def.F90",
    "lidort_lin_sup_brdf_def.F90",
    "lidort_lin_sup_ss_def.F90",
    "lidort_lin_sup_sleave_def.F90",
    "lidort_lin_sup_def.F90",
    "lidort_outputs_def.F90",
    "lidort_sup_brdf_def.F90",
    "lidort_sup_sleave_def.F90",
    "lidort_sup_ss_def.F90",
    "lidort_sup_def.F90",
    "lidort_inputs_def.F90",
        )

for lidort_def_fn in LIDORT_DEF_FILENAMES:
    DEF_FILENAMES.append(os.path.join(LIDORT_DEFS_DIR, lidort_def_fn))

# BRDF supplement drivers
BRDF_FILENAMES = (
    ("brdf_lin_sup_masters.F90", "brdf_linsup_masters_m"),
    ("brdf_sup_masters.F90", "brdf_sup_masters_m"),
    )

# Master routines
MASTER_FILENAMES = (
    "lidort_lcs_masters.F90",
    "lidort_lps_masters.F90",
    "lidort_inputs.F90",
    "lidort_masters.F90",
    )

# Preserve order for easier diffing against earlier versions of generation
MASTER_FILE_MODULES = OrderedDict() 

# BRDF masters do not follow the pattern of filename being the same as module name
for fn, module_name in BRDF_FILENAMES:
    MASTER_FILE_MODULES[os.path.join(BRDF_MAIN_DIR, fn)] = module_name

# Master modules have same module name as base filename
for fn in MASTER_FILENAMES:
    MASTER_FILE_MODULES[os.path.join(LIDORT_MAIN_DIR, fn)] = os.path.splitext(fn)[0]

# Add testing routines
MASTER_FILE_MODULES[os.path.join(LIDORT_BASE_DIR, 'lidort_test/lidort_sup_accessories.F90')] = "lidort_sup_accessories"

# Ignore these when parsing each file, they are loaded by the use statements
# into each file we parse
DEF_IGNORE_MODULES = ('lidort_pars') 

# Ignore these since they should be private inside LIDORT
MASTERS_IGNORE_ROUTINES = ( 'lidort_lcs_fourier_master',
                            'lidort_lps_fourier_master',
                            'lidort_fourier_master',
                            'lidort_init_inputs',
                            'lidort_read_inputs',
                            'lidort_check_input',
                            'lidort_check_input_thread',
                            'lidort_derive_input',
                            'lidort_sleave_input_checker', 
                          )

MASTERS_CONSTRUCTOR_ARGUMENTS={}
MASTERS_CONSTRUCTOR_ARGUMENTS["brdf_linsup_masters_m"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["brdf_sup_masters_m"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_lcs_masters"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_lps_masters"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_masters"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_input_master"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_sup_accessories"] = ( 'brdf_sup_in',
                                                            'lidort_fixin',
                                                            'lidort_modin', )

# These types will be passed instead of being
# made attributes
BRDF_TYPES = ( 'lidort_sup_inputs', 'lidort_linsup' )

NON_ATTRIBUTE_TYPES={}
NON_ATTRIBUTE_TYPES["lidort_lcs_master_module"] = BRDF_TYPES
NON_ATTRIBUTE_TYPES["lidort_lps_master_module"] = BRDF_TYPES
NON_ATTRIBUTE_TYPES["lidort_master_module"] = BRDF_TYPES

# Not sure if it really matters to supply this to fparser or not
USE_SOURCES = (
    PARS_FILENAME,
    )
PARS_IGNORE_VARIABLES = ( # Ignore debug write formats
                          "^dw", 
                          # Ignore kind storage variables
                          "lidort_dpkind",
                          "lidort_spkind",
                          "fpk",
                        )

##
# Output control

INTERFACE_TYPES_MODULE_NAME = "lidort_interface_types"
F_INTERFACE_TYPES_FILENAME = "%s.F90" % INTERFACE_TYPES_MODULE_NAME
F_INTERFACE_TYPES_IO_FILENAME = "%s_io.F90" % INTERFACE_TYPES_MODULE_NAME
F_INTERFACE_TYPES_IO_TMPL = "interface_types_io.f90"
H_INTERFACE_TYPES_FILENAME = "%s.h" % INTERFACE_TYPES_MODULE_NAME
I_INTERFACE_TYPES_FILENAME = "%s.i" % INTERFACE_TYPES_MODULE_NAME
TST_INTERFACE_TYPES_FILENAME = "%s_test.cc" % INTERFACE_TYPES_MODULE_NAME

INTERFACE_MASTERS_NAME = "lidort_interface_masters"
F_INTERFACE_MASTERS_FILENAME = "%s.F90" % INTERFACE_MASTERS_NAME
F_INTERFACE_MASTERS_IO_FILENAME = "%s_io.F90" % INTERFACE_MASTERS_NAME
F_INTERFACE_MASTERS_IO_TMPL = "interface_masters_io.f90"
H_INTERFACE_MASTERS_FILENAME = "%s.h" % INTERFACE_MASTERS_NAME
I_INTERFACE_MASTERS_FILENAME = "%s.i" % INTERFACE_MASTERS_NAME
TST_INTERFACE_MASTERS_FILENAME = "%s_test.cc" % INTERFACE_MASTERS_NAME

############################################################

def use_pars_var(var_name):
    # Do include certain pars variabls
    if not type(var_name) is str:
        return False
    elif re.match("fmt_", var_name):
        return False
    else:
        return True

# Use this for where else for fparse to search for information
source_only_list = [os.path.join(LIDORT_DEFS_DIR, fn) for fn in USE_SOURCES]

# Parse PARS file
pars_tree = fapi.parse(os.path.join(LIDORT_DEFS_DIR, PARS_FILENAME),
                       isfree=True, isstrict=False,
                       source_only=source_only_list)

pars_module = pars_tree.a.module[LIDORT_PARS_NAME]
pars_var_names = filter(use_pars_var, pars_module.a.variable_names)
pars_wrapper = TypeClass(pars_module, 
                         variable_names=pars_var_names, 
                         c_store_var_name="{name}",
                         ignore_variables=PARS_IGNORE_VARIABLES)

# Parse definition files
def_type_vars = []
for def_fn in DEF_FILENAMES:
    def_tree = fapi.parse(def_fn, isfree=True, isstrict=False,
                          source_only=source_only_list)
    for mod_name, mod_obj in def_tree.a.module.items():
        # Ignore some modules coming from use statements
        if mod_name not in DEF_IGNORE_MODULES:
            # Pull in the types from each module
            # Scan through content manually because .a.type_decls
            # does not preserve order and order is important 
            # when types are used within other types
            for content_obj in mod_obj.content:
                # Cannot link to private types
                if isinstance(content_obj, block.Type) and content_obj.is_public():
                    def_type_vars.append(content_obj)

# Specialized routines for copying supplement BRDF values into
# the Lidort BRDF types
def sup_brdf_copy_routines(type_obj):
    sup_config = {}
    if type_obj.name.find("linsup") >= 0:
        source_obj = "Brdf_Linsup_Outputs"
    else:
        source_obj = "Brdf_Sup_Outputs"

    copy_routine = """

void copy_from_sup({source_obj}& supp_obj) {{ 
  // This is only safe for LIDORT 3.6 because the BRDF and LIDORT
  // sup structures are the exact same structure, but with different
  // names. This MUST be reevaluated on future LIDORT versions
  void* sup_ptr = supp_obj.fortran_type_ptr();
  {name}_c_copy(&sup_ptr, &fortran_type_c);
}}
""".format(source_obj=source_obj, name=type_obj.name)

    sup_config['addl_c_public_methods'] = copy_routine

    return sup_config 

def_type_wrappers = []
for def_type_obj in def_type_vars:
    addl_config = {}
    if re.search('lidort_.*sup_brdf', def_type_obj.name):
        addl_config.update( sup_brdf_copy_routines(def_type_obj) )
    def_type_wrappers.append( TypeClass(def_type_obj, **addl_config) )

# Parse master files
def lidort_variable_configuration(routine_obj, var_obj):
    constructor_arg = False
    class_attribute = False

    if var_obj.name in MASTERS_CONSTRUCTOR_ARGUMENTS.get(routine_obj.parent.name, []):
        constructor_arg = True
        class_attribute = True
    elif type(var_obj.typedecl) is ftypes.Type and not var_obj.name in NON_ATTRIBUTE_TYPES.get(routine_obj.parent.name, []):
        class_attribute = True

    var_conf = { "constructor_arg": constructor_arg,
                 "class_attribute": class_attribute,
                 "store_bool_as_int": True,
               }

    if type(var_obj.typedecl) is ftypes.Type:
        var_conf["accessor_return_const"] = False
        var_conf["accessor_method_const"] = False 

    return var_conf

masters_list = parse_master_files(MASTER_FILE_MODULES,
                                  variable_configuration=lidort_variable_configuration,
                                  ignore_routines=MASTERS_IGNORE_ROUTINES)

# ############################################################

###
# Create C types interface 
write_cpp_type_classes(H_INTERFACE_TYPES_FILENAME, I_INTERFACE_TYPES_FILENAME, def_type_wrappers, pars_wrapper, base_class_name=DEF_BASE_CLASS_NAME, addl_includes=["<vector>"])

###
# Create C type tests
write_cpp_type_tests(TST_INTERFACE_TYPES_FILENAME, def_type_wrappers, pars_wrapper, suite_name=INTERFACE_TYPES_MODULE_NAME, addl_includes=[H_INTERFACE_TYPES_FILENAME])

###
# Create F type interface
write_fortran_type_module(F_INTERFACE_TYPES_FILENAME, def_type_wrappers, pars_wrapper, module_name=INTERFACE_TYPES_MODULE_NAME, addl_modules=[LIDORT_PARS_NAME])

write_fortran_type_module(F_INTERFACE_TYPES_IO_FILENAME, def_type_wrappers, pars_wrapper, module_name=INTERFACE_TYPES_MODULE_NAME, addl_modules=[LIDORT_PARS_NAME], template_fn=F_INTERFACE_TYPES_IO_TMPL)

###
# Create C masters classes
write_cpp_master_classes(H_INTERFACE_MASTERS_FILENAME, I_INTERFACE_MASTERS_FILENAME, 
                         masters_list,
                         h_template_fn="lidort_masters.h",
                         addl_includes=[H_INTERFACE_TYPES_FILENAME])

###
# Create F masters wrapper
write_fortran_master_modules(F_INTERFACE_MASTERS_FILENAME, 
                             masters_list, 
                             addl_modules=[LIDORT_PARS_NAME])

write_fortran_master_modules(F_INTERFACE_MASTERS_IO_FILENAME, 
                             masters_list, 
                             addl_modules=[LIDORT_INTERFACE_TYPES_IO_NAME],
                             template_fn=F_INTERFACE_MASTERS_IO_TMPL)
