// This file was auto-generated

%include "common.i"

%{
#include "{{ file_name }}.h"
%}

{%if pars_wrapper %}%fp_shared_ptr(FullPhysics::{{ pars_wrapper.class_name }});{% endif %}

namespace FullPhysics {

{% if pars_wrapper %}%nodefaultctor {{ pars_wrapper.class_name }};{% endif %}
%nodefaultctor {{ type_base_class_name }};

{% if pars_wrapper %}
struct {{ pars_wrapper.class_name }} {

  {% for var in pars_wrapper.variables -%}
  {{ var.c_store_var_decls(with_prefix=False, const=True)|c_decl_lines|indent(2) }}
  {% endfor %}
  static {{ pars_wrapper.class_name }}& instance();

};
{% endif %}

class {{ type_base_class_name }} {
public:
  void* fortran_type_ptr();

  std::string print_to_string() const;
};

{% for type in type_classes -%}
class {{ type.class_name }} : public {{ type_base_class_name }} {
public:
  {{ type.class_name }}();
  {{ type.class_name }}(const {{ type.class_name }}& src);
  ~{{ type.class_name }}();

  {% for var in type.variables -%}
  {{ var.c_get_accessor_sig() }};
  {% if not var.is_f_type("Character") -%}
  {{ var.c_set_accessor_sig() }};
  {% endif %}
  {% endfor %}

  {{ type.c_output_operator(use_accessor=True, as_printable=True).signature }};
};

{% endfor %}

}

