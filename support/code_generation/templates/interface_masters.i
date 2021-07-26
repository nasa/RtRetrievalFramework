// This file was auto-generated

%include "common.i"

{% from "constructor_def.macro" import constructor -%}

%{
#include "{{ file_name }}.h"
%}

{% if types_filename %}%import "{{ types_filename }}.i"{% endif %}

{% for master in master_classes -%}
%fp_shared_ptr(FullPhysics::{{master.class_name}});
{% endfor %}
namespace FullPhysics {

{% for master in master_classes %}
{% set def = master.c_class_def() %}
class {{ master.class_name }} {

public:
  {{ constructor(master.class_name, def.constructor_args) }};
  virtual ~{{ master.class_name }}();
  std::string print_to_string() const;

  {% for attr in master.attributes -%}
  {{ attr.i_get_accessor_sig()|indent(2) }}{% if "python_attribute" not in attr.i_get_accessor_sig() %};{% endif %}
  {% endfor -%}
  {% for routine in master.routines %}{%- set wrap = routine.c_wrapper() %}
  {{ wrap.signature }};{% endfor %}
};
{% endfor %}
}
