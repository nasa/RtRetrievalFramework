{% for master in master_classes %}
module {{ master.name|upper }}_WRAP

use iso_c_binding
use {{ master.name }}
{%- for mod_name in addl_used_modules %}
use {{ mod_name }}
{%- endfor %}
{%- for use_name in master.use_modules %}
use {{ use_name }}{% endfor %}

! This module was auto-generated 

implicit none

contains
{% for routine in master.routines %}
{% set wrap = routine.f_wrapper() -%}
{% set routine_beg = "subroutine " + wrap.wrapper_name -%}
{% set arg_indent = routine_beg|length + 2 -%}
{% set call_indent = wrap.wrapped_name|length + 8 -%}
{{ routine_beg }} ({{ wrap.argument_decls|f_arg_name_list|indent(arg_indent) }}) bind(C)
  {%- for use_name in routine.use_modules %}
  use {{ use_name }}{% endfor %}
  {%- if wrap.argument_decls|count > 0 %}

  ! Arguments
  {{ wrap.argument_decls|f_decl_lines()|indent(2) }}
  {%- endif %}
  {%- if wrap.local_decls|count > 0 %}

  ! Local variables
  {{ wrap.local_decls|f_decl_lines|indent(2) }}
  {%- endif %}
  {%- if wrap.input_conversion|count > 0 %}

  ! Convert input arguments
  {{ wrap.input_conversion|join("\n")|indent(2) }}
  {%- endif %}

  call {{ wrap.wrapped_name }}({{ wrap.call_arg_names|join(", &\n")|indent(call_indent) }})
  {%- if wrap.output_conversion|count > 0 %}

  ! Convert output arguments
  {{ wrap.output_conversion|join("\n")|indent(2) }}
  {%- endif %}

end subroutine {{ wrap.wrapper_name }}
{% endfor %}
end module {{ master.name|upper }}_WRAP
{% endfor %}
