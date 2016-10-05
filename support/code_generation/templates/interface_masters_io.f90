{% set lun_index = 500 -%}
{% for master in master_classes -%}
module {{ master.name|upper }}_IO

use iso_c_binding
{%- for mod_name in addl_used_modules %}
use {{ mod_name }}
{%- endfor %}
{%- for use_name in master.all_use_modules %}
use {{ use_name }}{% endfor %}

! This module was auto-generated 

implicit none

contains

{% for routine_type in ["read", "write"] -%}
{% set f_routine_name = master.name + "_" + routine_type -%} 
{% set routine_wrap = master.attributes_wrapper(f_routine_name, ignore="thread") -%}
{% set wrap = routine_wrap.f_wrapper() -%}

{% set routine_beg = "subroutine " + wrap.wrapper_name -%}
{% set arg_indent = routine_beg|length + 2 -%}
{% set call_indent = wrap.wrapped_name|length + 8 -%}
{{ routine_beg }} (filename_in, filename_in_len, {{ wrap.argument_decls|f_arg_name_list|indent(arg_indent) }}) bind(c)
  {%- if wrap.argument_decls|count > 1 %}

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  {{ wrap.argument_decls|f_decl_lines()|indent(2) }}
  {%- endif %}
  {%- if wrap.local_decls|count > 0 %}

  ! Local variables
  {{ wrap.local_decls|f_decl_lines|indent(2) }}
  {%- endif %}
  {%- if wrap.input_conversion|count > 0 %}
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  {{ wrap.input_conversion|join("\n")|indent(2) }}
  {%- endif %}
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call {{ wrap.wrapped_name }}(filename_lcl, {{ wrap.call_arg_names|join(", &\n")|indent(call_indent) }})
  {%- if wrap.output_conversion|count > 0 %}

  ! Convert output arguments
  {{ wrap.output_conversion|join("\n")|indent(2) }}
  {%- endif %}

end subroutine {{ wrap.wrapper_name }}

{% set routine_beg = "subroutine " + f_routine_name -%}
{% set arg_indent = routine_beg|length + 2 -%}
{{ routine_beg }} (filename, {{ wrap.argument_decls|f_arg_name_list|indent(arg_indent) }}) 
  ! Arguments
  character (len=*), intent(in) :: filename
  {% for attr in routine_wrap.arguments -%}
  {{ attr.f_arg_var_decls(use_type=True, intent="inout")|f_decl_lines|indent(2) }}
  {% endfor %}
  
  open ({{lun_index}}, file=filename, form="unformatted", access="sequential")
  {% for attr in routine_wrap.arguments -%}
  {% if attr.is_f_type("Type") -%}
  call {{ attr.typedecl.name }}_f_{{ routine_type }}({{lun_index}}, {{ attr.f_arg_var_name }})
  {% else -%}
  {{ routine_type }}(UNIT={{lun_index}}) attr.f_arg_var_name
  {% endif -%}
  {% endfor -%}
  close({{lun_index}})

end subroutine {{ f_routine_name }}

{% endfor %}
 
end module {{ master.name|upper }}_IO
{% set lun_index = lun_index + 1 %}
{% endfor %}
