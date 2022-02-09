{% from "constructor_def.macro" import constructor -%}
{% from "include_array.macro" import includes -%}
#ifndef {{ file_name|upper }}_H
#define {{ file_name|upper }}_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"
{{ includes(addl_includes) }}

/* This file was auto-generated */

namespace FullPhysics {

{% for master in master_classes -%}
{%- set def = master.c_class_def() -%}
//-----------------------------------------------------------------------
// Links to module: "{{ master.name }}" in file: "{{ master.source_basename }}"
//-----------------------------------------------------------------------

extern "C" {
{%- if has_read_write %}{%- for routine_type in ["read", "write"] %}
{%- set f_routine_name = master.name + "_" + routine_type %} 
{%- set routine_wrap = master.attributes_wrapper(f_routine_name, ignore="thread") %}
{%- set wrap = routine_wrap.c_wrapper() %}
  void {{ wrap.wrapper_name }}(const char* filename, const int* filename_len, {{ routine_wrap.extern_arg_decls()|join(", ") }});
{%- endfor %}{% endif %}
{%- for routine in master.routines %}
  void {{ routine.f_wrapper_name }}({{ routine.extern_arg_decls()|join(", ") }});
{%- endfor %}
}

class {{ master.class_name }} : public virtual GenericObject {

public:
  {{ constructor(master.class_name, def.constructor_args, def.class_inits) }} 
  { 
    {{ def.constructor_code|join("\n")|indent(4) }}
    // Initialize type pointers
    {% for attr in master.attributes if attr.is_f_type("Type") and not attr.constructor_arg -%}
    {{ attr.c_store_var_name }}.reset( new {{ attr.class_name() }}() );
    {% endfor %}
  }

  virtual ~{{ master.class_name }}() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  {% for attr in master.attributes -%}
  {% if not attr.accessor_return_const and not attr.accessor_method_const -%}
  {{ attr.c_get_accessor_code()|indent(2) }}
  {{ attr.c_get_accessor_code(return_const=True, method_const=True)|indent(2) }}
  {% else -%}
  {{ attr.c_get_accessor_code()|indent(2) }}
  {% endif -%}
  {% if attr.is_f_type("Type") -%}
  {{ attr.c_get_accessor_code(shared_ptr=True)|indent(2) }}
  {% endif -%}
  {% if not attr.is_f_type("Character") and attr.is_input_arg() and not attr.constructor_arg -%}
  {{ attr.c_set_accessor_code()|indent(2) }}
  {% endif %}

  {% endfor %}
  {% for routine in master.routines -%}
  {%- set wrap = routine.c_wrapper() -%}
  {{ wrap.signature }} {
    {% block pre_wrap scoped %}{{ wrap.pre_wrap_code|join("\n")|indent(4) }}{% endblock %}
    
    {{ wrap.f_routine_name }}({{ wrap.f_arguments|join(", ") }});
    {% block post_wrap scoped %}{{ wrap.post_wrap_code|join("\n")|indent(4) }}{% endblock %}
  }
{% endfor -%}
  {% if has_read_write %}{% for routine_type in ["read", "write"] -%}
  {% set f_routine_name = master.name + "_" + routine_type -%} 
  {% set routine_wrap = master.attributes_wrapper(f_routine_name, ignore="thread") -%}
  {% set wrap = routine_wrap.c_wrapper() %} 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void {{ routine_type }}_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    {{ wrap.pre_wrap_code|join("\n")|indent(4) }}

    {{ wrap.wrapper_name }}(filename_lcl, &filename_in_len, {{ wrap.f_arguments|join(", ") }});
    {{ wrap.post_wrap_code|join("\n")|indent(4) }}
  }
  {% endfor %}{% endif %}
  {%set out_op = master.c_output_operator(use_accessor=True) -%}
  {{ out_op.signature }} {
    {{ out_op.contents|indent(4) }}
  }

private:
  {{ def.attribute_decls|c_decl_lines|indent(2) }}
};

{% endfor %}

}
#endif
