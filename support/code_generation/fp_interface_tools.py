import re
import string
from string import Formatter

import collections

from fparser import typedecl_statements as ftypes, api as fapi

from jinja2 import Environment, FileSystemLoader

from variable_wrappers import *
from class_wrappers import *

################################################################################

def parse_master_files(filename_modules, **kwargs):
    master_objs = []
    for master_fn, mod_name in filename_modules.items():
        master_tree = fapi.parse(master_fn, isfree=True, isstrict=False)
        mod_obj = master_tree.a.module[mod_name]
        master_class_obj = MasterClass(mod_obj, **kwargs)
        master_objs.append(master_class_obj)
    return master_objs

################################################################################

def c_arg_decl_list(arg_dict):
    arg_list = []
    for name, decl in arg_dict.items():
        if isinstance(decl, collections.Iterable) and not isinstance(decl, str):
            var_str = "%s %s" % (decl[0], name)
        else:
            var_str = "%s %s" % (decl, name)
        arg_list.append(var_str)
    return ", ".join(arg_list)

def c_decl_lines(arg_dict):
    arg_lines = []
    for name, decl in arg_dict.items():
        if isinstance(decl, collections.Iterable) and not isinstance(decl, str):
            var_str = "%s %s%s;" % (decl[0], name, decl[1])
        else:
            var_str = "%s %s;" % (decl, name)
        arg_lines.append(var_str)
    return "\n".join(arg_lines)

def f_arg_name_list(arg_dict):
    return ", &\n".join([ "%s" % n for n in arg_dict.keys() ])

def f_decl_lines(arg_dict):
    arg_lines = []
    for name, decl in arg_dict.items():
        if isinstance(decl, collections.Iterable) and not isinstance(decl, str):
            var_str = "%s :: %s%s" % (decl[0], name, decl[1])
        else:
            var_str = "%s :: %s" % (decl, name)
        arg_lines.append(var_str)
    return "\n".join(arg_lines)

def one_or_more_lines(lines, postfix=""):
    if not isinstance(lines, collections.Iterable) or isinstance(lines, str):
        lines = [ lines ]
    return (postfix+"\n").join(lines)

TMPL_ENV = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))
TMPL_ENV.filters['c_arg_decl_list'] = c_arg_decl_list
TMPL_ENV.filters['c_decl_lines'] = c_decl_lines
TMPL_ENV.filters['f_arg_name_list'] = f_arg_name_list
TMPL_ENV.filters['f_decl_lines'] = f_decl_lines
TMPL_ENV.filters['one_or_more_lines'] = one_or_more_lines

def write_fortran_master_modules(filename, master_classes, template_fn="interface_masters.f90", addl_modules=[]):
    template = TMPL_ENV.get_template(template_fn)

    print("Writing to file: %s" % filename)
    with open(filename, 'w') as o_stream:
        o_stream.write(template.render(master_classes=master_classes,
                                       addl_used_modules=addl_modules))


def write_cpp_master_classes(h_filename, i_filename, master_classes, 
                             types_filename=None, has_read_write=True,
                             h_template_fn="interface_masters.h", 
                             i_template_fn="interface_masters.i",
                             addl_includes=[]):

    h_template = TMPL_ENV.get_template(h_template_fn)
    i_template = TMPL_ENV.get_template(i_template_fn)

    if types_filename is not None:
        full_types_filename = os.path.splitext(os.path.basename(types_filename))[0]
    else:
        full_types_filename = None

    for fn, tmpl in zip((h_filename, i_filename), (h_template, i_template)):
        print("Writing to file: %s" % fn)
        with open(fn, 'w') as o_stream:
            o_stream.write(tmpl.render(master_classes=master_classes,
                                       addl_includes=addl_includes,
                                       file_name=os.path.splitext(os.path.basename(fn))[0],
                                       types_filename=full_types_filename,
                                       has_read_write=has_read_write))

def write_cpp_type_tests(filename, type_classes, pars_wrapper=None, suite_name="interface_types", template_fn="interface_types_test.cc", addl_includes=[]):
    template = TMPL_ENV.get_template(template_fn)

    print("Writing to file: %s" % filename)
    with open(filename, 'w') as o_stream:
        o_stream.write(template.render(type_classes=type_classes,
                                       pars_wrapper=pars_wrapper,
                                       suite_name=suite_name,
                                       addl_includes=addl_includes))

def write_fortran_type_module(filename, type_classes, pars_wrapper=None, module_name="interface_types", template_fn="interface_types.f90", addl_modules=[]):
    template = TMPL_ENV.get_template(template_fn)

    print("Writing to file: %s" % filename)
    with open(filename, 'w') as o_stream:
        o_stream.write(template.render(type_classes=type_classes,
                                       pars_wrapper=pars_wrapper,
                                       module_name=module_name,
                                       addl_used_modules=addl_modules))

def write_cpp_type_classes(h_filename, i_filename, type_classes,
                           pars_wrapper=None,
                           h_template_fn="interface_types.h", 
                           i_template_fn="interface_types.i",
                           base_class_name="InterfaceStructure",
                           addl_includes=[]):

    h_template = TMPL_ENV.get_template(h_template_fn)
    i_template = TMPL_ENV.get_template(i_template_fn)

    for fn, tmpl in zip((h_filename, i_filename), (h_template, i_template)):
        print("Writing to file: %s" % fn)
        with open(fn, 'w') as o_stream:
            o_stream.write(tmpl.render(type_classes=type_classes,
                                       addl_includes=addl_includes,
                                       pars_wrapper=pars_wrapper,
                                       type_base_class_name=base_class_name,
                                       file_name=os.path.splitext(os.path.basename(fn))[0]))
