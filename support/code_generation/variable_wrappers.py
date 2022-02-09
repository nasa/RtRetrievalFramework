import re

from fparser import typedecl_statements as ftypes
from fparser.base_classes import Variable

from collections import OrderedDict, Iterable

##
# Type conversions

C_DECL_MAPPING = {
    ftypes.Integer: "{b_start}int{pointer}{b_end}",
    ftypes.Real: { 32: "{b_start}float{pointer}{b_end}", 
                   64: "{b_start}double{pointer}{b_end}" },
    ftypes.Character: "{b_start}char{pointer}{b_end}",
    ftypes.Logical: "{b_start}bool{pointer}{b_end}",
    ftypes.Type: "void*{pointer}",
    }

F2C_DECL_MAPPING = {
    ftypes.Integer: "integer(c_int){addl_decl}",
    ftypes.Real: { 32: "real(c_float){addl_decl}", 
                   64: "real(c_double){addl_decl}" },
    ftypes.Character: "character(kind=c_char, len={length}){addl_decl}",
    "argument_char": "character(kind=c_char) {addl_decl}",
    # Note that the logical should map to a int on the C side for pointers
    ftypes.Logical: { 8: "logical(c_bool){addl_decl}", 
                     32: "logical(kind=4){addl_decl}"},
    ftypes.Type: "type(c_ptr){addl_decl}",
    }

KIND_TYPE_BIT_MAPPING = { "lidort_spkind": 32,
                          "lidort_dpkind": 64,
                          "lidort_fpkind": 64,
                          "spk": 34,
                          "dpk": 64,
                          "fpk": 64,
                          "dp": 64,
                          "ffp": 64,
                          }

##
# Output control
FORTRAN_LINE_LIMIT = 132

# These names should not be used as bare names
RESERVED_NAMES = ( "alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand", "bitor", "bool",
                   "break", "case", "catch", "char", "char16_t", "char32_t", "class", "compl", "const",
                   "constexpr", "const_cast", "continue", "decltype", "default", "delete", "do", "double",
                   "dynamic_cast", "else", "enum", "explicit", "export", "extern", "false", "float",
                   "for", "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace",
                   "new", "noexcept", "not", "not_eq", "nullptr", "operator", "or", "or_eq", "private",
                   "protected", "public", "register", "reinterpret_cast", "return", "short", "signed",
                   "sizeof", "static", "static_assert", "static_cast", "struct", "switch", "template",
                   "this", "thread_local", "throw", "true", "try", "typedef", "typeid", "typename",
                   "union", "unsigned", "using(1)", "virtual", "void", "volatile", "wchar_t", "while",
                    "xor", "xor_eq", )

# Change any name that uses a RESERVED_NAME as such to avoid conflicts
RESERVED_NAME_CHANGE_TEMPLATE = "f_%s"

def get_decl_mapping(variable_obj, mapping_dict, **fmt_vals):
    if fmt_vals.get("srch_type", None) == None:
        srch_type = type(variable_obj.typedecl)
    else:
        srch_type = fmt_vals["srch_type"]
        
    try:
        mapped_val = mapping_dict[srch_type]
    except KeyError:
        raise KeyError("Could not find typedecl key object: %s among: %s" % (srch_type, mapping_dict.keys()))

    if type(mapped_val) is dict:
        var_kind = variable_obj.typedecl.get_kind()
        # Put declared types into KIND_TYPE_BIT_MAPPING
        # Which could probably be built up by further parsed
        # code analysis
        # Running get_bit_size() when the kinds is not an integer
        # in the source code will fail
        bit_size = fmt_vals.get("bit_size", None)
        if bit_size == None:
            bit_size = KIND_TYPE_BIT_MAPPING.get(var_kind, None)
        if bit_size == None:
            bit_size = variable_obj.typedecl.get_bit_size()
        mapped_val = mapped_val[bit_size]

    if "dim" not in fmt_vals:
        fmt_vals["fmt"] = ""
    fmt_vals["name"] = fmt_vals.get("name", None)
    if fmt_vals["name"] != None:
        raise DeprecationWarning("Name is deprecated, do not use")

    if variable_obj.init != None:
        fmt_vals["init"] = fix_var_init(variable_obj, variable_obj.init)

    fmt_vals["length"] = fmt_vals.get("length", None)
    if fmt_vals["length"] == None:
         fmt_vals["length"]= variable_obj.typedecl.get_length()

    mapped_val = mapped_val.format(**fmt_vals)

    return mapped_val

def fix_var_init(var_obj, init_str):

    if init_str is None:
        raise Exception("No initialization string for variable object:", var_obj)

    # Convert configure replacement value to preprocessor macro
    c_match = re.match("@([^@]+)@", init_str)
    if c_match:
        init_str = c_match.groups()[0].upper()

    # Change single quotes to doubles, since single quotes
    # in C/C++ means char not string
    init_str = re.sub("^[']", "\"", init_str)
    init_str = re.sub("[']$", "\"", init_str)

    # Remove _fpk type specifier from init strings
    init_str = re.sub("_fpk$", "", init_str)

    # Remove _ffp type specifier from init strings
    init_str = re.sub("_ffp$", "", init_str)

    # Remove _fpk type specifier from init strings
    init_str = re.sub("_dp$", "", init_str)

    # Change numbers with d0 into e0 
    init_str = re.sub("d0$", "e0", init_str)
    
    if type(var_obj.typedecl) is ftypes.Character:
        # Convert to initializer
        # Did have {0} here but that only works
        # with the C++0x standard. An empty
        # value passed as inside of a class
        # initalizer seems to be okay.
        init_str = re.sub('""', "", init_str)

    # Remove _fpk type specifier from init strings
    init_str = re.sub(".TRUE.", "true", init_str)
    init_str = re.sub(".FALSE.", "false", init_str)

    return init_str

def get_f_wrap_routine_name(module_name, routine_name, postfix="wrap", **kwargs):
    return_name = ""
    if module_name != None and len(module_name) > 0:
        return_name += module_name + "_"
    return_name += routine_name + "_" + postfix
    # Remove some extraneous strings from wrapper name
    return_name = return_name.replace("lidort_","").replace("module_","")
    return_name = return_name.replace("optimized_","")
    return return_name

def get_c_wrap_routine_name(module_name, routine_name):
    routine_name = routine_name.replace("lidort_","").replace("_plus","")
    routine_name = re.sub(".*inputmaster", "read_config", routine_name)
    routine_name = re.sub(".*input_master", "read_config", routine_name)
    routine_name = re.sub(".*master", "run", routine_name)
    return routine_name

def get_class_name(module_name):
    module_name = module_name.replace("_module","").title()
    module_name = re.sub("_[mM]$", "", module_name)
    module_name = re.sub("_optimized", "", module_name, flags=re.IGNORECASE)
    return module_name

class VariableWrapper(Variable):

    @classmethod
    def instance(cls, variable_obj, store_bool_as_int=False, char_in_type=False, **kwargs):
        """Dispatch the appropriate wrapper class"""

        if type(variable_obj.typedecl) is ftypes.Character:
            if char_in_type:
                return CharInTypeWrapper(variable_obj, **kwargs)
            elif kwargs.get("class_attribute", False):
                return CharAttributeWrapper(variable_obj, **kwargs)
            else:
                return CharWrapper(variable_obj, **kwargs)
        elif type(variable_obj.typedecl) is ftypes.Logical:
            if store_bool_as_int:
                return LogicalIntWrapper(variable_obj, **kwargs)
            else:
                return LogicalWrapper(variable_obj, **kwargs)
        elif type(variable_obj.typedecl) is ftypes.Type:
            return TypeWrapper(variable_obj, **kwargs)
        else:
            return VariableWrapper(variable_obj, **kwargs)

    def __init__(self, var_obj, **kwargs):
        self.__dict__ = var_obj.__dict__.copy()
           
        # Call init routines
        self.init_decls(**kwargs)

    def _init_name_attribute(self, var_name, default, **kwargs):
        # Set name variable attribute to value passed through
        # keywords as long as it isn't equal to None
        # If we do not find the value in keywords or
        # the passed value is None then use the default instead
        passed_val = kwargs.get(var_name, default)
        if passed_val == None:
            init_attr_value = default.format(**self.__dict__)
        else:
            init_attr_value = passed_val.format(**self.__dict__)

        if init_attr_value in RESERVED_NAMES:
            init_attr_value = RESERVED_NAME_CHANGE_TEMPLATE % init_attr_value

        setattr(self, var_name, init_attr_value) 

    def init_decls(self, **kwargs):
        self.routine_obj = kwargs.get("routine_obj", None)
        self.class_attribute = kwargs.get("class_attribute", False)
        self.constructor_arg = kwargs.get("constructor_arg", False)

        self.accessor_return_const = kwargs.get("accessor_return_const", True)
        self.accessor_method_const = kwargs.get("accessor_method_const", True)

        # Name of variable stored on the C side, generally as an attribute
        self._init_name_attribute("c_store_var_name", "{name}_", **kwargs)

        # Name of variable stored on the F side
        self._init_name_attribute("f_store_var_name", "{name}", **kwargs)

        # Name of variable when passing variable in C function arguments
        self._init_name_attribute("c_arg_var_name", "{name}_in", **kwargs)

        # Name of variable when passing variable in F function arguments
        self._init_name_attribute("f_arg_var_name", "{name}_in", **kwargs)

        # Name of variable converted to another type for sending to fortran
        self._init_name_attribute("c_conv_var_name", "{name}_lcl", **kwargs)

        # Name of variable converted inside of fortran wrapper
        self._init_name_attribute("f_conv_var_name", "{name}_lcl", **kwargs)

        # Container prefix for variables
        self._init_name_attribute("c_container_prefix", "", **kwargs)
        self._init_name_attribute("f_container_prefix", "", **kwargs)

        # Name of variable converted inside of fortran wrapper
        self._init_name_attribute("accessor_name", "{name}", **kwargs)


    def f_arg_var_decls(self, **kwargs):
        decls = OrderedDict()
        decls[self.f_arg_var_name] = self.f_declaration(add_intent=True, **kwargs)
        return decls

    def f_store_var_decls(self, with_prefix=True, **kwargs):
        decls = OrderedDict()
        if with_prefix:
            store_name = self.f_container_prefix + self.f_store_var_name
        else:
            store_name = self.f_store_var_name
        decls[store_name] = self.f_declaration(**kwargs)
        return decls

    def f_supporting_var_decls(self):
        supp_vars = OrderedDict()
        if self.is_array():
            supp_vars["%s_f_shapes" % self.name] = "integer(c_int), dimension(%d)" % (len(self.get_array_spec()))

        if self.is_f_type(ftypes.Character):
            supp_vars["%s_f_len" % self.name] = "integer(c_int)"
        else:
            supp_vars["%s_f_byte_size" % self.name] = "integer(c_int)"
        return supp_vars

    def f_var_comment(self):
        if self.is_array():
            return "dimension(" + self.f_dimension_str() + ")"
        else:
            return "scalar"

    def f_local_var_decls(self):
        return OrderedDict()

    def f_call_wrapper_names(self):
        return self.f_arg_var_decls().keys()

    def f_copy_from_c_code(self):
        return []

    def f_copy_to_c_code(self):
        return []

    def f_wrapper_routine_name(self):
        return get_f_wrap_routine_name(self.parent.parent.name, self.name, "get")

    def c_arg_var_decls(self, pointer=False, const=True, **kwargs):
        decls = OrderedDict()
        if self.is_array() and not pointer:
            decls[self.c_arg_var_name] = self.blitz_type(use_bool=True, const=const, **kwargs) + "&"
        else:
            decls[self.c_arg_var_name] = self.c_declaration(pointer=pointer, const=const, **kwargs) + "&"
        return decls

    def c_store_var_decls(self, pointer=False, with_prefix=True, **kwargs):
        decls = OrderedDict()
        if not self.is_array() and with_prefix:
            store_name = self.c_container_prefix + self.c_store_var_name
        else:
            store_name = self.c_store_var_name
        if self.is_array() and not pointer:
            decls[store_name] = self.blitz_type(**kwargs)
        else:
            decls[store_name] = self.c_declaration(pointer=pointer, **kwargs)
        return decls

    def c_extern_arg_decls(self, const=True, **kwargs):
        decls = OrderedDict()
        decls[self.c_arg_var_name] = self.c_declaration(const=const, pointer=True, bool_override=True, **kwargs)
        return decls

    def c_call_f_wrapper_args(self):
        if self.class_attribute:
            call_name = self.c_store_var_name
        else:
            call_name = self.c_arg_var_name    

        if self.is_array():
            return [ "%s.dataFirst()" % call_name ]
        else:
            return [ "&%s" % call_name ]

    def c_constructor_code(self):
        def get_storage_name(var_name):
            if var_name.isdigit():
                return var_name
            elif self.routine_obj != None and hasattr(self.routine_obj, "size_variables"):
                var_obj = self.routine_obj.size_variables.get(var_name, None)
                if var_obj != None:
                    return var_obj.c_store_var_name
                else:
                    return var_name
            else:
                return var_name

        init_code = []
        if self.is_array() and self.class_attribute:
            size_args = []
            for spec in self.get_array_spec():
                if len(spec) > 1:
                    size_range = get_storage_name(spec[1]) + "-" + get_storage_name(spec[0]) + "+1"
                    size_args.append(size_range)
                else:
                    size_args.append(get_storage_name(spec[0]))

            init_code += [ 
                "{var_name}.reference( {blitz_type}({size_args}, blitz::ColumnMajorArray<{dim}>()) );". \
                    format(var_name=self.c_store_var_name, 
                           blitz_type=self.blitz_type(),
                           size_args=", ".join(size_args),
                           dim=self.num_dim()) ]

        # Don't set something passed through constructor back to zero!
        if not self.constructor_arg:
            init_code.append("{var_name} = {zero_value};".format(var_name=self.c_store_var_name, zero_value=self.c_zero_value()))
                           
        return init_code

    def c_accessor_set_decl(self, **kwargs):
        return list(self.c_arg_var_decls(**kwargs).items())[0]

    def c_accessor_return_decl(self, const=False, **kwargs):
        const_str = "const " if const else ""
        (return_name, return_type) = list(self.c_store_var_decls(**kwargs).items())[0]
        return (return_name, const_str + return_type)

    def c_routine(self, signature, contents_arr, indent=2):
        routine_code = ""
        if not isinstance(signature, Iterable) or isinstance(signature, str):
            signature = [ signature ]
        for curr_sig in signature:
            routine_code += "{signature} {{\n{contents}\n}}\n".\
                            format(signature=curr_sig,
                                   contents="\n".join([(" "*indent) + cline for cline in contents_arr]))
        return routine_code

    def c_pre_wrap_code(self):
        return []

    def c_post_wrap_code(self):
        return []

    def c_get_accessor_sig(self, accessor_name=None, return_ref=True, return_const=None, method_const=None, **kwargs):
        # Default value if not overridden
        if return_const == None:
            return_const = self.accessor_return_const 
        if method_const == None:
            method_const = self.accessor_method_const

        return_name, return_type = self.c_accessor_return_decl(const=return_const, **kwargs)
        return "{out_var_type}{ref_str} {accessor_name}(){const}". \
            format(out_var_type=return_type,
                   ref_str = "&" if return_ref else "",
                   accessor_name=accessor_name or self.accessor_name,
                   const = " const" if method_const else "")

    def i_get_accessor_sig(self, accessor_name=None, return_ref=True, return_const=None, **kwargs):

        return_name, return_type = self.c_accessor_return_decl(const=return_const, **kwargs)
        return "%python_attribute({accessor_name}, {out_var_type}{ref_str})". \
            format(out_var_type=return_type,
                   ref_str = "&" if return_ref else "",
                   accessor_name=accessor_name or self.accessor_name)

    def c_set_accessor_sig(self, accessor_name=None, set_const=True, **kwargs):
        set_name, set_type = self.c_accessor_set_decl(const=set_const, **kwargs)
        return "void {accessor_name}({in_var_type} {set_var_name})". \
            format(accessor_name=accessor_name or self.accessor_name,
                   in_var_type=set_type, 
                   set_var_name=set_name)

    def c_get_accessor_code(self, accessor_name=None, **kwargs):
        return_name, return_type = self.c_accessor_return_decl(**kwargs)

        get_code = [ "return %s;" % return_name ]

        return self.c_routine(self.c_get_accessor_sig(accessor_name, **kwargs), get_code)

    def c_set_accessor_code(self, accessor_name=None, **kwargs):
        arg_decls = self.c_arg_var_decls()
        store_decls = self.c_store_var_decls()

        set_code = []
        for a_name, s_name in zip(arg_decls.keys(), store_decls.keys()):
            set_code.append("%s = %s;" % (s_name, a_name))
    
        return self.c_routine(self.c_set_accessor_sig(accessor_name), set_code)

    def f_dimension_str(self, no_extent=False):
        if no_extent:
            dim_decl_str = ",".join(":" * len(self.get_array_spec()))
        else:
            dim_strs = []
            for dim_decl in self.get_array_spec():
                dim_strs.append( ":".join(dim_decl).upper() )
            dim_decl_str = ", ".join(dim_strs)

        return dim_decl_str

    def c_output_operator_code(self, disp_var_name=None, container_name="", use_accessor=False, pre_var_output = "", post_var_output = ""):
        if disp_var_name == None:
            disp_var_name = self.name

        if use_accessor:
            get_call = self.accessor_name + "()"
        else:
            get_call = self.c_store_var_name
        
        if self.is_array():
            pre_var_output += " << std::endl"

        return """\n  << "{disp_var_name}: "{pre_var_output} << {container_name}{get_call} {post_var_output} << std::endl""".format(disp_var_name=disp_var_name, get_call=get_call, pre_var_output=pre_var_output, post_var_output=post_var_output, container_name=container_name), False
                       
    def f_intent_str(self, intent=None):
        if self.intent or intent:
            intent_str = "intent(%s)" % (intent if intent else self.intent[0].lower())
            return intent_str
        else:
            return None

    def f_declaration(self, add_dimension=True, no_dim_extent=False, f_pointer=False, c_pointer=False, add_intent=False, **fmt_vals):
        fmt_vals["addl_decl"] = ""

        dim_str = ""
        if self.is_array() and (add_dimension or c_pointer):
            dim_str = "dimension(" + self.f_dimension_str(no_extent=no_dim_extent) + ")"
            if not c_pointer:
                fmt_vals["addl_decl"] += ", " + dim_str

        if f_pointer:
            fmt_vals["addl_decl"] += ", pointer"
        elif c_pointer:
            fmt_vals["srch_type"] = ftypes.Type
           
        intent_str = self.f_intent_str(fmt_vals.get("intent", None))
        if add_intent and intent_str != None:
            fmt_vals["addl_decl"] += ", " + intent_str 

        f_decl = get_decl_mapping(self, F2C_DECL_MAPPING, **fmt_vals)
    
        # Wrap at first space before line limit
        if len(f_decl) > FORTRAN_LINE_LIMIT:
            space_loc = f_decl.rfind(" ", 0, FORTRAN_LINE_LIMIT-2)
            f_decl = f_decl[0:space_loc] + "&\n" + f_decl[space_loc:]
    
        return f_decl.strip()

    def c_supporting_vars(self):
        supp_vars = []
        if self.is_array():
            supp_vars.append("int %s_f_shapes[%d];" % (self.c_store_var_name, len(self.get_array_spec())))

        if self.is_f_type(ftypes.Character):
            supp_vars.append("int %s_f_len;" % self.c_store_var_name)
        else:
            supp_vars.append("int %s_f_byte_size;" % self.c_store_var_name)

        return supp_vars

    def c_declaration(self, const=False, pointer=False, blitz_array=False, bool_override=False, **fmt_vals):
        # For pointers to logicals, use integer on the C side
        if type(self.typedecl) is ftypes.Logical and (pointer or blitz_array) and not bool_override:
            fmt_vals["srch_type"] = ftypes.Integer

        if pointer:
            fmt_vals["pointer"] = "*"
        else:
            fmt_vals["pointer"] = ""

        fmt_vals["b_start"] = ""
        fmt_vals["b_end"] = ""
        if blitz_array and self.is_array():
            blitz_fmt = {}
            self.blitz_type(blitz_fmt=blitz_fmt)
            fmt_vals["b_start"] = blitz_fmt["b_start"]
            fmt_vals["b_end"] = blitz_fmt["b_end"]

        fmt_vals["type_name"] = fmt_vals.get("type_name", None)
        if fmt_vals["type_name"] == None:
            fmt_vals["type_name"] = self.class_name()

        c_decl = get_decl_mapping(self, C_DECL_MAPPING, **fmt_vals).strip()

        if const:
            c_decl = "const " + c_decl

        return c_decl

    def blitz_type(self, const=False, use_bool=True, num_dims=None, blitz_fmt=None, **kwargs):
        if blitz_fmt == None:
            blitz_fmt = {}

        if num_dims == None:
            num_dims = len(self.get_array_spec())

        blitz_fmt["b_start"] = "blitz::Array<"
        blitz_fmt["b_end"] = ", %d>" % num_dims
        blitz_fmt["pointer"] = ""
        blitz_fmt["comment"] = ""

        # Use integer on the C side for fortran logicals
        if type(self.typedecl) is ftypes.Logical and not use_bool:
            blitz_fmt["srch_type"] = ftypes.Integer

        map_str = get_decl_mapping(self, C_DECL_MAPPING, **blitz_fmt)

        if const:
            map_str = "const " + map_str

        return map_str

    def is_f_type(self, type_spec):
        if isinstance(type_spec, str):
            type_spec = getattr(ftypes, type_spec)
        return type(self.typedecl) is type_spec

    def intent_str(self):
        if self.intent != None and len(self.intent) > 0:
            return self.intent[0].lower()
        else:
            return "inout"

    def is_input_arg(self):
        return re.search("in", self.intent_str())

    def is_output_arg(self):
        return re.search("out", self.intent_str())

    def init_string(self, object_name, other_vars=[]):
        var_names = [ v.name for v in other_vars ]
        init_str = fix_var_init(self, self.init)
        self_refs = re.findall("([a-z][a-z_1-9]+)", init_str)
        if self_refs != None:
            for ref_name in self_refs:
                if ref_name in var_names:
                    init_str = re.sub("("+ref_name+")","%s.\g<0>" % object_name, init_str)
        return init_str

    def expected_dims(self):
        expt_dims = []
        for curr_dim in self.get_array_spec():
            dim_val = ""
            if len(curr_dim) == 2:
                minus_str = ""
                if curr_dim[0] != "0":
                    minus_str = "-%s" % curr_dim[0]
                dim_val = "%s%s+1" % (curr_dim[1], minus_str)
            else:
                dim_val = "%s" % curr_dim[0]
            expt_dims.append(dim_val)
        return expt_dims

    def f_zero_value(self):
        # Not all types have the ability to query the zero value:
        typedecl_by_name = self.typedecl.get_type_decl(self.typedecl.name)
        if isinstance(self.typedecl, ftypes.Type) and typedecl_by_name is None:
            raise Exception("Could not determine f_zero_value for: %s in routine %s" % (self.name, self.routine_obj.name))
        else:
            return self.typedecl.get_zero_value()

    def c_zero_value(self):
        # Not all types have the ability to query the zero value:
        typedecl_by_name = self.typedecl.get_type_decl(self.typedecl.name)
        if isinstance(self.typedecl, ftypes.Type) and typedecl_by_name is None:
            if self.routine_obj is not None:
                raise Exception("Could not determine c_zero_value for: %s in routine %s" % (self.name, self.routine_obj.name))
        else:
            return fix_var_init(self, self.typedecl.get_zero_value())

    def f_lower_bounds(self):
         if self.is_array():
             lower_bounds = "(" + ",".join(["&\n    lbound(fortran_type_f%%%s,%d)"%(self.f_store_var_name,idx+1) for idx in range(len(self.get_array_spec()))]) + ")"
         else:
             lower_bounds = ""
         return lower_bounds

    def c_array_len_str(self):
        # Use char definition if [] if * used in fortran
        if self.typedecl.get_length() == "*":
            length_str = "[]"
        else:
            length_str = "[%s]" % self.typedecl.get_length()
        return length_str


    def num_dim(self):
        if self.is_array():
            return len(self.get_array_spec())
        else:
            return 0

    def class_name(self):
        return get_class_name(self.typedecl.name) 
 
############################################################

class CharWrapper(VariableWrapper):
    def __init__(self, var_obj, **kwargs):

        if type(var_obj.typedecl) is not ftypes.Character:
            raise Exception("Wrapped object must of type: %s" % ftypes.Character)

        VariableWrapper.__init__(self, var_obj, **kwargs)

        self._init_name_attribute("f_shape_var_name", "{f_arg_var_name}_shape_{{0}}", **kwargs)
        # In case the passed f_shape_var is just the size of the array and
        # we doe not want to waste time copying maybe uninitialized values
        self._init_name_attribute("f_shape_loop_var_name", "{f_shape_var_name}", **kwargs)
        self._init_name_attribute("c_shape_var_name", "{c_arg_var_name}_shape_{{0}}", **kwargs)
        self._init_name_attribute("f_len_var_name", "{f_arg_var_name}_len", **kwargs)
        self._init_name_attribute("c_len_var_name", "{c_arg_var_name}_len", **kwargs)

        # For local variables in copying code
        self._init_name_attribute("f_dim_index_var_name", "dim_idx_{{0}}", **kwargs)
        self._init_name_attribute("f_len_index_var_name", "len_idx", **kwargs)
        self._init_name_attribute("f_lower_bound_var_name", "lb_{{0}}", **kwargs)

        # Since we do not return chars directly
        self._init_name_attribute("c_return_var_name", "{name}_ret", **kwargs)

    def f_arg_var_decls(self, **kwargs):
        decls = OrderedDict()

        # Define variables storing the shape of array
        shape_name_list = []
        if self.is_array():
            for dim_idx in range(len(self.get_array_spec())):
                dim_var_name = self.f_shape_var_name.format(dim_idx+1)
                shape_name_list.append(dim_var_name)
                decls[dim_var_name] = "integer(c_int), intent(in)"

        # Define length variable
        decls[self.f_len_var_name] = "integer(c_int), intent(in)"

        # Define variable containing character data itself
        dim_str = ""
        if self.is_array():
            dim_str = ", ".join(shape_name_list) + ", "
        length_str = dim_str + self.f_len_var_name + "+1"
        decls[self.f_arg_var_name] = (self.f_declaration(add_intent=True, intent="inout", add_dimension=False, srch_type="argument_char"), "(%s)" % length_str)

        return decls

    def f_local_var_decls(self, conv_var=True):
        decls = OrderedDict()

        # Define a declaration to use for local variable for copying
        # to make it available for other code in the routine being created
        if conv_var:
            decls[self.f_conv_var_name] = self.f_declaration(length=self.f_len_var_name)

        for array_dim in range(1,self.num_dim()+1):
            decls[self.f_dim_index_var_name.format(array_dim)] = "integer"
            decls[self.f_lower_bound_var_name.format(array_dim)] = "integer"

        # All character variables have an associated length
        decls[self.f_len_index_var_name] = "integer"

        return decls

    def f_call_wrapper_names(self):
        return [ self.f_conv_var_name ]

    def f_copy_loop_start(self, f_var_name):
        loop_start = []
        # Use a variable for use as a loop index for the dimension
        # and to get the lower bound of the fortran array to
        # ensure we do not assume anything about how that
        # variable's extents are set up
        for array_dim in range(self.num_dim()):
            # Set up code for each dimension loop
            lower_bound_code =  \
                ("  "*(array_dim-1)) + \
                "{lb_var} = lbound({container}{var},{dim})".\
                format(lb_var=self.f_lower_bound_var_name.format(array_dim+1), 
                       dim=array_dim+1, var=f_var_name, 
                       container=self.f_container_prefix)
            loop_start.append(lower_bound_code)

            dim_loop_start = \
                ("  "*(array_dim-1)) + \
                "do {idx_var} = 1, {shape_var}". \
                format(idx_var=self.f_dim_index_var_name.format(array_dim+1), 
                       shape_var=self.f_shape_loop_var_name.format(array_dim+1))
        
            loop_start.append(dim_loop_start)

        # Start of length loop
        len_loop_start = ("  "*self.num_dim()) + \
            "do %s = 1, %s" % (self.f_len_index_var_name, self.f_len_var_name)
        loop_start.append(len_loop_start)

        return loop_start

    def f_copy_loop_end(self):
        # Set up list for end loop code so we can append after the fact
        # and append to from the inside out of the loop structure
        end_loop_code = []
        end_loop_code.append(("  "*self.num_dim())+"end do")
        
        for dim_idx in range(self.num_dim()):
            end_loop_code.append(("  "*dim_idx)+"end do")

        return end_loop_code

    def f_copy_loop_index_strings(self):
        # Set up indexing string for accessing data
        f_idx_str = "("
        c_idx_str = "("
        if self.num_dim() > 0:
            f_idx_str += ", ".join([ "%s-1+%s" % (self.f_dim_index_var_name.format(idx+1), self.f_lower_bound_var_name.format(idx+1)) for idx in range(self.num_dim())]) + ")("
            c_idx_str += ", ".join([ self.f_dim_index_var_name.format(idx+1) for idx in range(self.num_dim())]) + ", "

        all_idx_str = f_idx_str + ":)"
        f_idx_str += "%s:%s)" % (self.f_len_index_var_name, self.f_len_index_var_name)
        c_idx_str += "%s)" % self.f_len_index_var_name

        return (f_idx_str, c_idx_str, all_idx_str)

    def f_copy_to_c_code(self, source_name=None, dest_name=None, bounds_name=None):

        if bounds_name is None:
            if self.is_output_arg(): 
                bounds_name = self.f_conv_var_name
            else:
                bounds_name = self.f_arg_store_name

        copy_to_code = self.f_copy_loop_start(bounds_name)
        end_loop_code = self.f_copy_loop_end()

        if source_name == None:
            source_name = self.f_conv_var_name
        if dest_name == None:
            dest_name = self.f_arg_var_name

        f_idx_str, c_idx_str, all_idx_str = self.f_copy_loop_index_strings()

        copy_to_code.append(
            "{indent}{c_name}{c_idx_str} = &\n{indent}{indent}{f_name}{f_idx_str}". \
                format(c_name=dest_name, f_name=source_name, 
                       c_idx_str=c_idx_str, f_idx_str=f_idx_str, indent="  "*(self.num_dim()+1)) )

        # Append inner loop ending statement
        copy_to_code.append(end_loop_code[0])

        # Add null character
        copy_to_code.append(
            "{indent}len_idx = len_trim({f_name}{all_idx_str})+1". \
                format(f_name=source_name, all_idx_str=all_idx_str, 
                       indent="  "*self.num_dim(), container=self.f_container_prefix) )
        copy_to_code.append(
            "{indent}{c_name}{c_idx_str} = c_null_char". \
                format(c_name=dest_name, c_idx_str=c_idx_str, indent="  "*self.num_dim()) )
  
        # Append outer loop ending statements
        copy_to_code += end_loop_code[1:]

        return copy_to_code

    def f_copy_from_c_code(self, source_name=None, dest_name=None):

        bounds_name = self.f_arg_var_name

        copy_from_code = self.f_copy_loop_start(bounds_name)
        end_loop_code = self.f_copy_loop_end()

        if source_name == None:
            source_name = self.f_arg_var_name
        if dest_name == None:
            dest_name = self.f_conv_var_name

        f_idx_str, c_idx_str, all_idx_str = self.f_copy_loop_index_strings()

        # Define code for copying of values from one array to another
        copy_from_code.append( 
            "{indent}{f_name}{f_idx_str} = &\n{indent}{indent}{c_name}{c_idx_str}". \
                format(c_name=source_name, f_name=dest_name, 
                       c_idx_str=c_idx_str, f_idx_str=f_idx_str, 
                       indent="  "*(self.num_dim()+1), container=self.f_container_prefix) )

        # Append inner loop ending statement
        copy_from_code.append(end_loop_code[0])
  
        # Append outer loop ending statements
        copy_from_code += end_loop_code[1:]

        return copy_from_code

    def c_interface_type(self):
        """Type exposed to C++ for use in accessors or through a C++ function call"""

        # Build type based on number of dimensions of string
        interface_type = "std::string"
        for array_dim in range(self.num_dim()):
            interface_type = "std::vector< %s >" % interface_type
        
        return interface_type

    def c_extern_arg_decls(self, **kwargs):
        decls = OrderedDict()

        # Define variables describing shape of array
        for dim_idx in range(self.num_dim()):
            dim_var_name = self.c_shape_var_name.format(dim_idx+1)
            decls[dim_var_name] = self.c_declaration(const=True, pointer=True, srch_type=ftypes.Integer, **kwargs) 

        # Define length variable
        decls[self.c_len_var_name] = self.c_declaration(const=True, pointer=True, srch_type=ftypes.Integer)

        # Define data variable itself
        decls[self.c_arg_var_name] = self.c_declaration(const=True, pointer=True, length="")
        return decls

    def c_arg_var_decls(self, pointer=False, **kwargs):
        decls = OrderedDict()
        decls[self.c_arg_var_name] = "const " + self.c_interface_type() + "&"
        return decls

    def c_accessor_set_decl(self, **kwargs):
        return (self.c_arg_var_name, "const " + self.c_interface_type() + "&")

    def c_get_accessor_sig(self, accessor_name=None):
        return VariableWrapper.c_get_accessor_sig(self, accessor_name, return_ref=False)

    def i_get_accessor_sig(self, accessor_name=None):
        return self.c_get_accessor_sig(accessor_name=accessor_name)

    def c_call_f_wrapper_args(self):
        call_args = []

        for dim_idx in range(self.num_dim()):
            call_args.append(self.c_shape_var_name.format(dim_idx+1))

        call_args += [ "&%s" % self.c_len_var_name,
                       "%s" % self.c_conv_var_name ]

        return call_args

    def c_pre_wrap_code(self):
        set_code = []
        # Setting simples strings to pass directly, arrays of vectors not supported
        if not self.is_array():
            set_code += [
                "const char* {ret_var} = {set_var}.c_str();".\
                    format(ret_var=self.c_conv_var_name, set_var=self.c_arg_var_name),
                "int {len_var} = (int) {set_var}.size();".\
                    format(len_var=self.c_len_var_name, set_var=self.c_arg_var_name), ]
        return set_code

    def c_copy_from_char_array(self, src_name, dest_name):
        copy_code = []

        # Note.. this will not handle more than 1 dimension
        assignment_str = " = "
        arr_index_str = ""
        for dim_idx in range(self.num_dim()):
            dim_index_var = "dim_%d_idx" % dim_idx
            assignment_str = ".push_back"
            arr_index_str += "%s, " % dim_index_var
            copy_code.append(
                "  "*dim_idx + \
                    "for(int {idx_var} = 0; {idx_var} < {lcl_var}.extent({dim_idx}); {idx_var}++)".\
                    format(idx_var=dim_index_var, lcl_var=src_name, 
                           dim_idx=dim_idx))

        copy_code.append(
            "  "*self.num_dim() + \
                "{ret_var}{assignment}( std::string(std::string({conv_var}({arr_index}blitz::Range::all()).begin(), {conv_var}({arr_index}blitz::Range::all()).end()).c_str()) );". \
                format(ret_var=dest_name, assignment=assignment_str, 
                       conv_var=src_name, arr_index=arr_index_str))

        return copy_code

    def c_output_operator_code(self, disp_var_name=None, container_name="", use_accessor=False, pre_var_output="", post_var_output=""):
        if use_accessor:
            get_call = self.accessor_name + "()"
        else:
            get_call = self.c_store_var_name

        if self.is_array():
            return """
  << "{disp_var_name}: " << std::endl;
std::vector< std::string > {lcl_var_name} = {container_name}{get_call};
for(unsigned int idx = 0; idx < {lcl_var_name}.size(); idx++)
  if ( {lcl_var_name}[idx].length() > 0 )
    output_stream << "  [" << idx << "]: \\"" << {lcl_var_name}[idx] << "\\"" << std::endl;""".format(disp_var_name=disp_var_name, get_call=get_call, pre_var_output=pre_var_output, container_name=container_name, lcl_var_name=self.c_conv_var_name), True

        else:
            pre_var_output += ' << "\\""'
            post_var_output += '<< "\\""'
            return VariableWrapper.c_output_operator_code(self, disp_var_name, container_name, use_accessor, pre_var_output, post_var_output)

class CharAttributeWrapper(CharWrapper):
    """When wrapping a character type stored as a class attribute
    where the space allocation happens on the C++ side"""

    def num_chars(self, dim_idx):
        return self.typedecl.selector[dim_idx]

    def char_len(self):
        if self.typedecl.get_length() == '*':
            raise Exception("Can not get character length of char(len=*) for variable named: %s in routine: %s" % (self.name, self.routine_obj.name))
        return int(self.typedecl.get_length())+1

    def c_constructor_code(self):
        # Size our stored blitz array as the shape variables plus
        # the character length
        size_spec = [ self.num_chars(d) for d in range(self.num_dim()) ] + [ self.char_len() ]

        init_code = [ 
            "{stor_var}.reference( {stor_type}({size_spec}, blitz::ColumnMajorArray<{num_dim}>()) );".\
                format(stor_var=self.c_store_var_name,
                       stor_type=self.c_store_var_decls()[self.c_store_var_name],
                       size_spec=", ".join([str(s) for s in size_spec]),
                       num_dim=self.num_dim()+1),
            "{stor_var} = '\\0';".format(stor_var=self.c_store_var_name),
            ]


        return init_code

    def c_store_var_decls(self, **kwargs):
        decls = VariableWrapper.c_store_var_decls(self, **kwargs)
        if self.is_array():
            decls[self.c_store_var_name] = self.blitz_type(num_dims=self.num_dim()+1)
        else:
            decls[self.c_store_var_name] = self.blitz_type(num_dims=1)
        return decls


    def c_pre_wrap_code(self):
        set_code = []
        # Setting simples strings to pass directly, arrays of vectors not supported
        for dim_idx in range(self.num_dim()):
            set_code += [ 
                "int {shape_var} = (int) {set_var}.extent({shape_idx});".\
                    format(shape_var=self.c_shape_var_name.format(dim_idx+1), 
                           set_var=self.c_store_var_name,
                           shape_idx=dim_idx) ]

        # Don't count the null character when passing size to Fortran,
        # it will account for it internally
        set_code += [
            "int {len_var} = (int) {set_var}.extent({len_idx}) - 1;".\
                format(len_var=self.c_len_var_name, 
                       set_var=self.c_store_var_name,
                       len_idx=self.num_dim()) ]

        return set_code

    def c_call_f_wrapper_args(self):
        call_args = []

        for dim_idx in range(self.num_dim()):
            call_args.append("&%s" % self.c_shape_var_name.format(dim_idx+1))

        call_args += [ "&%s" % self.c_len_var_name,
                       "%s.dataFirst()" % self.c_store_var_name ]

        return call_args

    def c_accessor_return_decl(self, const=False, **kwargs):
        const_str = "const " if const else ""
        return (self.c_return_var_name, self.c_interface_type())

    def c_set_accessor_code(self, accessor_name=None):
        raise NotImplementedError("Setting of character attributes not yet supported.")
    
    def c_get_accessor_code(self, accessor_name=None, **kwargs):
        return_name, return_type = self.c_accessor_return_decl()
        get_code = []
        get_code.append("%s %s;" % (return_type, return_name))

        get_code += self.c_copy_from_char_array(self.c_store_var_name, self.c_return_var_name)
            
        get_code.append("return %s;" % return_name)

        return self.c_routine(self.c_get_accessor_sig(accessor_name, **kwargs), get_code)


class CharInTypeWrapper(CharWrapper):
    """When wrapping a character type stored as in a user defined type
    Where the allocation happens on the fotran side"""

    def c_store_var_decls(self, **kwargs):
        decls = VariableWrapper.c_store_var_decls(self, **kwargs)
        if self.is_array():
            decls[self.c_store_var_name] = decls.get(self.c_store_var_name, None)
        else:
            decls[self.c_store_var_name] = (decls.get(self.c_store_var_name, None), self.c_array_len_str())
        return decls

    def c_accessor_return_decl(self, const=False, **kwargs):
        const_str = "const " if const else ""
        return (self.c_return_var_name, const_str + self.c_interface_type())

    def c_get_accessor_code(self, accessor_name=None, **kwargs):
        return_name, return_type = self.c_accessor_return_decl(**kwargs)

        conv_var_name = self.c_conv_var_name
        conv_var_type = self.blitz_type(num_dims=self.num_dim()+1)

        get_code = []

        get_code.append("%s %s;" % (return_type, return_name))

        size_arguments = []
        if self.is_array():
            size_arguments += ["%s%s_f_shapes[%d]" % (self.c_container_prefix, self.c_store_var_name, dim_idx) for dim_idx in range(self.num_dim())]
        size_arguments.append("%s%s_f_len" % (self.c_container_prefix, self.c_store_var_name))

        size_references = ["&%s" % arg for arg in size_arguments]

        get_code.append(
            "{conv_var_type} {conv_var_name} = {conv_var_type}({size_arguments}+1, blitz::ColumnMajorArray<{num_blitz_dim}>());". \
                format(conv_var_type=conv_var_type,
                       conv_var_name=conv_var_name,
                       size_arguments=", ".join(size_arguments),
                       num_blitz_dim=len(size_arguments)) )

        get_code.append(
            "{wrapper_routine_name}(const_cast<void**>(&fortran_type_c), {size_references}, {conv_var_name}.dataFirst());". \
                format(wrapper_routine_name=self.f_wrapper_routine_name(),
                       size_references=", ".join(size_references),
                       conv_var_name=conv_var_name) )

        get_code += self.c_copy_from_char_array(conv_var_name, self.c_return_var_name)
            
        get_code.append("return %s;" % return_name)

        return self.c_routine(self.c_get_accessor_sig(accessor_name, **kwargs), get_code)

    def c_set_accessor_code(self, accessor_name=None):
        set_code = self.c_pre_wrap_code()
        return self.c_routine(self.c_set_accessor_sig(accessor_name), set_code)


############################################################

class LogicalWrapper(VariableWrapper):

    def __init__(self, var_obj, **kwargs):

        if type(var_obj.typedecl) is not ftypes.Logical:
            raise Exception("Wrapped object must of type: %s" % ftypes.Logical)

        VariableWrapper.__init__(self, var_obj, **kwargs)

    def f_arg_var_decls(self, **kwargs):
        decls = OrderedDict()
        decls[self.f_arg_var_name] = self.f_declaration(add_intent=True, bit_size=8, **kwargs)
        return decls

    def f_local_var_decls(self):
        decls = OrderedDict()
        decls[self.f_conv_var_name] = self.f_declaration()
        return decls

    def f_call_wrapper_names(self):
        return [ self.f_conv_var_name ]

    def f_copy_from_c_code(self):
        return ["{f_name} = {c_name}".format(f_name=self.f_conv_var_name, 
                                             c_name=self.f_arg_var_name)]

    def f_copy_to_c_code(self):
        return ["{c_name} = {f_name}".format(f_name=self.c_conv_var_name, 
                                             c_name=self.f_arg_var_name)]
    def f_zero_value(self):
        return ".FALSE."

    def c_zero_value(self):
        return "false"

class LogicalIntWrapper(LogicalWrapper):

    def c_arg_var_decls(self, const=True, **kwargs):
        return VariableWrapper.c_arg_var_decls(self, const=const, **kwargs)

    def c_extern_arg_decls(self, const=True, **kwargs):
        return VariableWrapper.c_extern_arg_decls(self, const=const, **kwargs)

    def c_get_accessor_sig(self, accessor_name=None, **kwargs):
        # Do not allow return by reference due to need to make a copy of
        # fortran values
        return VariableWrapper.c_get_accessor_sig(self, accessor_name=None, return_ref=False, **kwargs)

    def c_get_accessor_code(self, accessor_name=None, **kwargs):
        if self.is_array():
            get_code = [ 
                "blitz::Array<bool,{num_dim}> as_bool({var_name}.shape());". \
                    format(var_name=self.c_store_var_name, 
                           num_dim=self.num_dim()),
                "as_bool = blitz::where({var_name} != 0, true, false);". \
                    format(var_name=self.c_store_var_name),
                "return as_bool;" ]
        else:
            return_name, return_type = self.c_accessor_return_decl()
            get_code = [ "return %s != 0;" % return_name ]

        return self.c_routine(self.c_get_accessor_sig(accessor_name, **kwargs), get_code)


    def c_set_accessor_code(self, accessor_name=None):
        if self.is_array():
            set_code = [ 
                "blitz::Array<int,{num_dim}> as_int({var_name}.shape());". \
                    format(var_name=self.c_store_var_name, 
                           num_dim=self.num_dim()),
                "as_int = blitz::where({set_var_name} == true, FORTRAN_TRUE_INT, 0);".\
                    format(set_var_name=self.c_arg_var_name),
                "{var_name} = as_int;".format(var_name=self.c_store_var_name) ]
        else:
            arg_decls = self.c_arg_var_decls()
            store_decls = self.c_store_var_decls()

            set_code = []
            for a_name, s_name in zip(arg_decls.keys(), store_decls.keys()):
                set_code.append("%s = %s ? FORTRAN_TRUE_INT : 0;" % (s_name, a_name))

        return self.c_routine(self.c_set_accessor_sig(accessor_name), set_code)


############################################################

class TypeWrapper(VariableWrapper):
    def __init__(self, var_obj, **kwargs):
        
        if type(var_obj.typedecl) is not ftypes.Type:
            raise Exception("Wrapped object must of type: %s" % ftypes.Type)

        VariableWrapper.__init__(self, var_obj, **kwargs)

        self.shared_ptr_tmpl = "boost::shared_ptr<%s>"
        self.type_ptr_tmpl = "type({type_name}){intent}, pointer"

    def c_constructor_code(self):
        # No default initialization
        return []

    def f_arg_var_decls(self, use_type=False, **kwargs):
        decls = OrderedDict()
        if use_type:
            intent_str = ", " + self.f_intent_str(**kwargs)
            decls[self.f_arg_var_name] = self.type_ptr_tmpl.format(type_name=self.typedecl.name,
                    intent=intent_str)
        else:
            decls[self.f_arg_var_name] = self.f_declaration(add_intent=True, **kwargs)
        return decls

    def f_local_var_decls(self):
        decls = OrderedDict()
        decls[self.f_conv_var_name] = \
            self.type_ptr_tmpl.format(type_name=self.typedecl.name, intent="")
        return decls

    def f_call_wrapper_names(self):
        return [ self.f_conv_var_name ]

    def f_copy_from_c_code(self):
        return ["call c_f_pointer({c_name}, {f_name})".\
                    format(c_name=self.f_arg_var_name, f_name=self.f_conv_var_name)]

    def c_arg_var_decls(self, with_prefix=True, shared_ptr=False, **kwargs):
        decls = OrderedDict()
        if shared_ptr:
            decls[self.c_arg_var_name] = (self.shared_ptr_tmpl + "&") % self.class_name()
        else:
            decls[self.c_arg_var_name] = self.class_name() + "&"
        return decls

    def c_extern_arg_decls(self, const=False, **kwargs):
        return VariableWrapper.c_extern_arg_decls(self, const=const, **kwargs)

    def c_store_var_decls(self, with_prefix=True, **kwargs):
        decls = OrderedDict()
        decls[self.c_store_var_name] = self.shared_ptr_tmpl % self.class_name()
        return decls

    def c_accessor_set_decl(self, **kwargs):
        return (self.c_arg_var_name, self.class_name() + "&")

    def c_accessor_return_decl(self, const=False, shared_ptr=False, **kwargs):
        const_str = "const " if const else ""
        if shared_ptr:
            return_decl = const_str + self.shared_ptr_tmpl % self.class_name()
        else:
            return_decl = const_str + self.class_name()
        return (self.c_store_var_name, return_decl)

    def c_get_accessor_sig(self, accessor_name=None, shared_ptr=False, **kwargs):
        if shared_ptr:
            accessor_name = (accessor_name if accessor_name else self.accessor_name) + "_ptr"
        return VariableWrapper.c_get_accessor_sig(self, accessor_name, shared_ptr=shared_ptr, **kwargs) 

    def c_get_accessor_code(self, accessor_name=None, shared_ptr=False, **kwargs):
        return_name, return_type = self.c_accessor_return_decl(shared_ptr=shared_ptr)

        if shared_ptr:
            get_code = [ "return %s;" % return_name ]
        else:
            get_code = [ "return *%s;" % return_name ]

        return self.c_routine(self.c_get_accessor_sig(accessor_name, shared_ptr=shared_ptr, **kwargs), get_code)

    def c_set_accessor_code(self, accessor_name=None):
        arg_decls = self.c_arg_var_decls()
        store_decls = self.c_store_var_decls()

        set_code = []
        for a_name, s_name in zip(arg_decls.keys(), store_decls.keys()):
            set_code.append("void* src_ptr = %s.fortran_type_ptr();" % a_name)
            set_code.append("void* dst_ptr = %s->fortran_type_ptr();" % s_name)
            set_code.append("{name}_c_copy(&src_ptr, &dst_ptr);".format(name=self.typedecl.name))
    
        return self.c_routine(self.c_set_accessor_sig(accessor_name), set_code)


    def c_pre_wrap_code(self):
        # The fortran_type_ptr methods should defined for the Type wrapping class
        # so we can convert it a void pointer
        if self.class_attribute:
            set_var_name = self.c_store_var_name
        else:
            set_var_name = self.c_arg_var_name

        return [ "void* {container}{local_var_name} = {set_var_name}->fortran_type_ptr();".\
                     format(local_var_name=self.c_conv_var_name, 
                            container=self.c_container_prefix, 
                            set_var_name=set_var_name) ]

    def c_call_f_wrapper_args(self):
        return [ "&%s" % self.c_conv_var_name ]
