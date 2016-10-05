{% from "constructor_def.macro" import constructor -%}
// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This file was auto-generated

%{
#include "{{ file_name }}.h"
%}

namespace FullPhysics {

{% for master in master_classes %}
{% set def = master.c_class_def() %}
class {{ master.class_name }} {

public:
  {{ constructor(master.class_name, def.constructor_args) }};
  
  {% for attr in master.attributes -%}
  {{ attr.c_get_accessor_sig()|one_or_more_lines(";")|indent(2) }};
  {% endfor -%}
  {% for routine in master.routines %}{%- set wrap = routine.c_wrapper() %}
  {{ wrap.signature }};{% endfor %}
};
{% endfor %}
}
