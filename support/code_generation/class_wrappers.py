import os
import re

from io import StringIO
from collections import OrderedDict

from fparser import typedecl_statements as ftypes
from fparser import statements 

from variable_wrappers import *

###########################################################################

class WrapperClass(object):
    def __init__(self, wrap_obj=None, use_modules=None, **kwargs):
        self.class_name = ""
        self.wrap_obj = wrap_obj

        # Add any use statements from module or it's parent
        if use_modules:
            self.use_modules = use_modules
        else:
            self.use_modules = []

        if wrap_obj != None:
            for content_obj in wrap_obj.content + wrap_obj.parent.content:
                if isinstance(content_obj, statements.Use):
                    self.use_modules.append(content_obj.name)

        # Copy these kwargs as attributes, no processing needed
        for copy_name in ('addl_c_private_methods', 'addl_c_public_methods', 'addl_f_routines'):
            value = kwargs.get(copy_name, None)
            if value:
                setattr(self, copy_name, kwargs[copy_name])
            else:
                setattr(self, copy_name, "")

    def c_output_operator(self, variables, use_accessor=False, as_printable=False):
        decl_details = OrderedDict()
        max_var_len = max([ len(var.name) for var in variables ])

        if as_printable:
            decl_details["signature"] = "virtual void print(std::ostream &output_stream) const"
        else:
            decl_details["signature"] = "friend std::ostream& operator<<(std::ostream &output_stream, const {class_name} &obj)".format(class_name=self.class_name)

        decl_details["contents"] = ""
        decl_details["contents"] += """output_stream << "{class_name}:" << std::endl""".format(class_name=self.class_name)

        restart_stream = False
        for var_obj in variables:
            disp_var_name = " " * (max_var_len - len(var_obj.name)) + var_obj.name

            obj_name = ""
            if not as_printable:
                obj_name = "obj."

            if restart_stream:
                decl_details["contents"] += "\noutput_stream"

            var_contents, restart_stream = var_obj.c_output_operator_code(disp_var_name, obj_name, use_accessor)

            decl_details["contents"] += var_contents

        if not restart_stream:
            decl_details["contents"] += ";\n"
        if not as_printable:
            decl_details["contents"] += "return output_stream;\n"
    
        return decl_details;

######
       
class TypeClass(WrapperClass):
    def __init__(self, def_type_obj, variable_names=None, ignore_variables=None, f_container_name=None, store_orig_name=False, **kwargs):
        WrapperClass.__init__(self, def_type_obj, **kwargs)

        self.variables = []

        # Use all variables if specific variable list is not specified
        if variable_names == None:
            variable_names = def_type_obj.a.variable_names

        # Filter out unwanted variable names
        if ignore_variables:
            filtered_var_names = []
            for in_var_name in variable_names:
                ignored = False
                for ignore_re in ignore_variables:
                    if re.search(ignore_re, in_var_name):
                        ignored = True

                if not ignored:
                    filtered_var_names.append(in_var_name)

            variable_names = filtered_var_names

        if f_container_name == None:
            self.f_container_name = "fortran_type_f"

        for var_name in variable_names:
            var_obj = def_type_obj.get_variable(var_name)

            if type(var_obj.typedecl) is ftypes.Character:
                c_container_prefix="transfer_struct_c."
            else:
                c_container_prefix="*transfer_struct_c."

            wrap_obj = VariableWrapper.instance(var_obj, 
                                                c_container_prefix=c_container_prefix,
                                                f_container_prefix=self.f_container_name + "%",
                                                store_bool_as_int=True,
                                                char_in_type=True,
                                                **kwargs)
            self.variables.append(wrap_obj)

        self.f_type_name = def_type_obj.name
        self.f_parent_name = def_type_obj.parent.name
        self.f_bind_name = "%s_c" % self.f_type_name

        self.class_name = get_class_name(self.f_type_name)

        self.source_filename = def_type_obj.reader.file.name
        self.source_basename = os.path.basename(self.source_filename)

        self.class_inits = [ "%s(%s)" % (var.name, var.c_zero_value()) for var in self.variables ]

    def c_output_operator(self, use_accessor=False, as_printable=False):
        return WrapperClass.c_output_operator(self, self.variables, use_accessor, as_printable)

######

class RoutineWrapper(WrapperClass):
                                       
    def __init__(self, arguments, name, module_name=None, parent=None, variable_configuration=None, c_routine_name_func=get_c_wrap_routine_name, f_routine_name_func=get_f_wrap_routine_name, **kwargs):
        self.name = name
        self.module_name = module_name
        self.parent = parent
        self.c_routine_name_func = c_routine_name_func
        self.f_routine_name_func = f_routine_name_func
        self.f_wrapper_name = self.f_routine_name_func(module_name, self.name, **kwargs)

        self.size_variables = OrderedDict()
        self.arguments = []
        for var_obj in arguments: 

            var_conf = OrderedDict()
            if variable_configuration:
                if isinstance(variable_configuration, dict):
                    var_conf = variable_configuration
                else:
                    var_conf.update( variable_configuration(self, var_obj) )
            
            if isinstance(var_obj, VariableWrapper):
                wrap_obj = var_obj
            else:
                wrap_obj = VariableWrapper.instance(var_obj, routine_obj=self, **var_conf)
            self.arguments.append(wrap_obj)

            if var_conf.get("size_variable", False):
                self.size_variables[wrap_obj.name] = wrap_obj

    def extern_arg_decls(self):
        arg_decls = []
        for var_obj in self.arguments:
            for name, decl in var_obj.c_extern_arg_decls().items():
                arg_decls.append("%s %s" % (decl, name))
        return arg_decls
 
    def c_wrapper(self):
        c_args_decl = []
        f_args_list = []
        pre_wrap_code = []
        post_wrap_code = []
        for wrap_obj in self.arguments:
            # If the variable is an attribute then it should be
            # get/set and converted through accessor functions
            if not wrap_obj.class_attribute:
                for name, decl in wrap_obj.c_arg_var_decls().items():
                    c_args_decl.append("%s %s" % (decl, name))

            pre_wrap_code += wrap_obj.c_pre_wrap_code()
            post_wrap_code += wrap_obj.c_post_wrap_code()

            f_args_list += wrap_obj.c_call_f_wrapper_args()


        c_routine_name = self.c_routine_name_func(self.module_name, self.name)
        f_routine_name = self.f_routine_name_func(self.module_name, self.name)
        wrap_signature = "void {routine_name}({arguments})".format(routine_name=c_routine_name, arguments=", ".join(c_args_decl))

        return { "wrapper_name": self.f_wrapper_name,
                 "signature": wrap_signature,
                 "c_routine_name": c_routine_name,
                 "f_routine_name": f_routine_name,
                 "f_arguments": f_args_list,
                 "pre_wrap_code": pre_wrap_code,
                 "post_wrap_code": post_wrap_code,
                 }

    def f_wrapper(self):
        """Convenience routine to consolidate argument and local declarations"""
    
        # Output input types
        argument_decls = OrderedDict()
        local_decls = OrderedDict()
        pre_convert_code = []
        call_arg_names = []
        post_convert_code = []
        for wrap_obj in self.arguments:
            for name, decl in wrap_obj.f_arg_var_decls().items():
                argument_decls[name] = decl

            for name, decl in wrap_obj.f_local_var_decls().items():
                local_decls[name] = decl

            # Derived types will be passed as pointers            
            if wrap_obj.is_input_arg() or wrap_obj.is_f_type(ftypes.Type):
                pre_convert_code += wrap_obj.f_copy_from_c_code()

            call_arg_names += wrap_obj.f_call_wrapper_names()

            if wrap_obj.is_output_arg():
                post_convert_code += wrap_obj.f_copy_to_c_code()

        return { "wrapper_name": self.f_wrapper_name,
                 "wrapped_name": self.name,
                 "argument_decls": argument_decls,
                 "local_decls": local_decls,
                 "input_conversion": pre_convert_code,
                 "call_arg_names": call_arg_names,
                 "output_conversion": post_convert_code 
                 }
 
####################################################################################
    
class MasterClass(WrapperClass):

    def __init__(self, module_obj, ignore_routines=None, **kwargs):
        WrapperClass.__init__(self, module_obj, **kwargs)

        self.routines = []
        self.all_use_modules = []
           
        self.name = module_obj.name
        self.class_name = get_class_name(self.name)
        self.source_filename = module_obj.reader.file.name
        self.source_basename = os.path.basename(self.source_filename)

        for mod_routine in sorted(module_obj.a.module_provides.values(), key=lambda r: r.name):
            if mod_routine.is_public() and not mod_routine.name in ignore_routines:
                routine = Master_Routine(mod_routine, parent=self, **kwargs)
                self.routines.append(routine)
                self.all_use_modules += routine.use_modules

        self.attributes = []
        self.constructor_args = []

        seen_attrs = []
        seen_const = []
        for routine_obj in self.routines:
            for arg in routine_obj.arguments:
                if arg.class_attribute and not arg.name in seen_attrs:
                    self.attributes.append(arg)
                    seen_attrs.append(arg.name)
                if arg.constructor_arg and not arg.name in seen_const:
                    self.constructor_args.append(arg)
                    seen_const.append(arg.name)

    def attributes_wrapper(self, name, ignore=None, **kwargs):
        used_attrs = []
        if ignore != None:
            if type(ignore) is str:
                ignore = [ ignore ]
            for attr in self.attributes:
                for ig_pattern in ignore:
                    if not re.search(ig_pattern, attr.name):
                        used_attrs.append(attr)
        else:
            used_attrs = self.attributes

        return RoutineWrapper(used_attrs, name, **kwargs)

    def c_class_def(self):
    
        # Define constructor
        constr_arg_list = []
        init_list = []
        for var_obj in self.constructor_args:
            arg_decls = var_obj.c_arg_var_decls(shared_ptr=True)
            store_decls = var_obj.c_store_var_decls()

            for (a_name, a_decl), (s_name, s_decl) in zip(arg_decls.items(), store_decls.items()):
                constr_arg_list.append("%s %s" % (a_decl, a_name))
                init_list.append("%s(%s)" % (s_name, a_name))

        attribute_decls = OrderedDict()
        constructor_code = []
        for var_obj in self.attributes:
            for (s_name, s_decl) in var_obj.c_store_var_decls().items():
                attribute_decls[s_name] = s_decl
            constructor_code += var_obj.c_constructor_code()

        return { "constructor_args": constr_arg_list,
                 "constructor_code": constructor_code,
                 "class_inits": init_list,
                 "attribute_decls": attribute_decls }

    def c_output_operator(self, use_accessor=False, as_printable=False):
        return WrapperClass.c_output_operator(self, self.attributes, use_accessor, as_printable)

######
       
class Master_Routine(RoutineWrapper):
                                       
    def __init__(self, routine_obj, parent=None, variable_configuration=None, **kwargs):
        WrapperClass.__init__(self, routine_obj, **kwargs)

        arguments = [ routine_obj.get_variable(argname) for argname in routine_obj.args ]
        RoutineWrapper.__init__(self, arguments, routine_obj.name, routine_obj.parent.name, parent,
                variable_configuration, **kwargs)

