#!/usr/bin/env python

import os

# Requires F2PY G3 code from http://code.google.com/p/f2py/
# Used changeset:   81:69b8613d30cb
from fparser import api as fapi

from fp_interface_tools import *

try:
    # Python 2.7+
    from collections import OrderedDict
except ImportError:
    # http://pypi.python.org/pypi/ordereddict
    from ordereddict import OrderedDict

##
# Fortran input filenames

LIDORT_DEFS_DIR = os.path.join(os.environ["L2_EXE_PATH"], "thirdparty/lidort-3.5/lidort_def")
LIDORT_MAIN_DIR = os.path.join(os.environ["L2_EXE_PATH"], "thirdparty/lidort-3.5/lidort_main")

PARS_FILENAME = "lidort_pars.F90.in"
LIDORT_PARS_NAME = "lidort_pars"

DEF_BASE_CLASS_NAME = "Lidort_Structure"

# Derived type definitions
DEF_FILENAMES = (
    "brdf_sup_inputs_def.F90", 
    "lidort_brdf_sup_def.F90", 
    "lidort_inputs_def.F90", 
    "lidort_lc_outputs_def.F90", 
    "lidort_lin_inputs_def.F90", 
    "lidort_lp_outputs_def.F90", 
    "lidort_ls_brdf_sup_def.F90", 
    "lidort_ls_outputs_def.F90", 
    "lidort_outputs_def.F90",
    )

# Master routines
MASTER_FILENAMES = (
    "lidort_brdf_master_module.F90",
    "lidort_lcs_master_module.F90",
    "lidort_lps_master_module.F90",
    "lidort_ls_brdf_master_module.F90",
    "lidort_master_module.F90",
    )

# Preserve order for easier diffing against earlier versions of generation
MASTER_FILE_MODULES = OrderedDict() 
for fn in MASTER_FILENAMES:
    MASTER_FILE_MODULES[os.path.join(LIDORT_MAIN_DIR, fn)] = os.path.splitext(fn)[0]

# Ignore these when parsing each file, they are loaded by the use statements
# into each file we parse
DEF_IGNORE_MODULES = ('lidort_pars', 'lidort_typekinds')

# Ignore these since they should be private inside LIDORT
MASTERS_IGNORE_ROUTINES = ( 'brdf_quadrature_gaussian',
                            'brdf_quadrature_trapezoid',
                            'brdf_maker',
                            'brdf_function',
                            'brdf_fourier',
                            'brdf_maker_plus',
                            'brdf_function_plus',
                            'brdf_ls_fourier',
                            'lidort_lcs_fourier_master',
                            'lidort_lps_fourier_master',
                            'lidort_fourier_master')

MASTERS_CONSTRUCTOR_ARGUMENTS={}
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_lcs_master_module"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_lps_master_module"] = ( 'thread' )
MASTERS_CONSTRUCTOR_ARGUMENTS["lidort_master_module"] = ( 'thread' )

# These types will be passed instead of being
# made attributes
BRDF_TYPES = ( 'lidort_brdf_inputs', 'lidort_lsbrdf_inputs' )

NON_ATTRIBUTE_TYPES={}
NON_ATTRIBUTE_TYPES["lidort_lcs_master_module"] = BRDF_TYPES
NON_ATTRIBUTE_TYPES["lidort_lps_master_module"] = BRDF_TYPES
NON_ATTRIBUTE_TYPES["lidort_master_module"] = BRDF_TYPES

# Not sure if it really matters to supply this to fparser or not
USE_SOURCES = (
    PARS_FILENAME,
    "lidort_type_kinds.F90",
    )


##
# Output control

INTERFACE_TYPES_MODULE_NAME = "lidort_interface_types"
F_INTERFACE_TYPES_FILENAME = "%s.F90" % INTERFACE_TYPES_MODULE_NAME
H_INTERFACE_TYPES_FILENAME = "%s.h" % INTERFACE_TYPES_MODULE_NAME
I_INTERFACE_TYPES_FILENAME = "%s.i" % INTERFACE_TYPES_MODULE_NAME
TST_INTERFACE_TYPES_FILENAME = "%s_test.cc" % INTERFACE_TYPES_MODULE_NAME

INTERFACE_MASTERS_NAME = "lidort_interface_masters"
F_INTERFACE_MASTERS_FILENAME = "%s.F90" % INTERFACE_MASTERS_NAME
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
                         c_store_var_name="{name}")

# Parse definition files
def_type_vars = []
for def_fn in DEF_FILENAMES:
    def_tree = fapi.parse(os.path.join(LIDORT_DEFS_DIR, def_fn),
                                 isfree=True, isstrict=False,
                                 source_only=source_only_list)
    for mod_name, mod_obj in def_tree.a.module.items():
        # Ignore some modules coming from use statements
        if mod_name not in DEF_IGNORE_MODULES:
            # Pull in the types from each module
            for def_type_name, def_type_obj in mod_obj.a.type_decls.items():
                # Cannot link to private types
                if def_type_obj.is_public():
                    def_type_vars.append(TypeClass(def_type_obj))

# Parse master files
def lidort_variable_configuration(routine_obj, var_obj):
    constructor_arg = False
    class_attribute = False

    if var_obj.name in MASTERS_CONSTRUCTOR_ARGUMENTS.get(routine_obj.parent.name, []):
        constructor_arg = True
        class_attribute = True
    elif type(var_obj.typedecl) is ftypes.Type and not var_obj.name in NON_ATTRIBUTE_TYPES.get(routine_obj.parent.name, []):
        class_attribute = True

    return { "constructor_arg": constructor_arg,
             "class_attribute": class_attribute,
             "store_bool_as_int": True,
             }

masters_list = parse_master_files(MASTER_FILE_MODULES,
                                  variable_configuration=lidort_variable_configuration,
                                  ignore_routines=MASTERS_IGNORE_ROUTINES)

# ############################################################

###
# Create C types interface 
write_cpp_type_classes(H_INTERFACE_TYPES_FILENAME, I_INTERFACE_TYPES_FILENAME, def_type_vars, pars_wrapper, base_class_name=DEF_BASE_CLASS_NAME, addl_includes=["<vector>"])

###
# Create C type tests
write_cpp_type_tests(TST_INTERFACE_TYPES_FILENAME, def_type_vars, pars_wrapper, suite_name=INTERFACE_TYPES_MODULE_NAME, addl_includes=[H_INTERFACE_TYPES_FILENAME])

###
# Create F type interface
write_fortran_type_module(F_INTERFACE_TYPES_FILENAME, def_type_vars, pars_wrapper, module_name=INTERFACE_TYPES_MODULE_NAME, addl_modules=[LIDORT_PARS_NAME])

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
