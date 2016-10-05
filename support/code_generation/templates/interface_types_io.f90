module {{ module_name }}_io

use iso_c_binding
{%- for mod_name in addl_used_modules %}
use {{ mod_name }}
{%- endfor %}

! This module was auto-generated 

implicit none

contains
{% for type in type_classes %}
! Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
! Allocs and initializes type
subroutine {{ type.f_type_name }}_c_write(lun, fortran_type_c) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type({{ type.f_type_name }}), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call {{ type.f_type_name }}_f_write(lun, fortran_type_f)

end subroutine {{ type.f_type_name }}_c_write

subroutine {{ type.f_type_name }}_f_write(lun, fortran_type_f) 
  use {{ type.f_parent_name }}
  {% for module in type.use_modules -%}
  use {{ module }}
  {% endfor %}
  integer, intent(in) :: lun
  type({{ type.f_type_name }}), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  {% for var in type.variables if var.is_f_type("Type") -%}
  {{ var.f_local_var_decls()|f_decl_lines|indent(2) }}  
  {% endfor %}
  ! Get pointer to types
  {% for var in type.variables if var.is_f_type("Type") -%}
  {{ var.f_conv_var_name }} => fortran_type_f%{{ var.name }}
  {% endfor %}
  {% for var in type.variables -%}
  {% if var.is_f_type("Type") -%}
  call {{ var.typedecl.name }}_f_write(lun, {{ var.f_conv_var_name }})
  {% else -%}
  write(UNIT=lun) fortran_type_f%{{ var.name }}
  {% endif -%}
  {% endfor %}
end subroutine {{ type.f_type_name }}_f_write

subroutine {{ type.f_type_name }}_c_read(lun, fortran_type_c) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type({{ type.f_type_name }}), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call {{ type.f_type_name }}_f_read(lun, fortran_type_f)

end subroutine {{ type.f_type_name }}_c_read

subroutine {{ type.f_type_name }}_f_read(lun, fortran_type_f) 
  use {{ type.f_parent_name }}
  {% for module in type.use_modules -%}
  use {{ module }}
  {% endfor %}
  integer, intent(in) :: lun
  type({{ type.f_type_name }}), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  {% for var in type.variables if var.is_f_type("Type") -%}
  {{ var.f_local_var_decls()|f_decl_lines|indent(2) }}  
  {% endfor %}
  ! Get pointer to types
  {% for var in type.variables if var.is_f_type("Type") -%}
  {{ var.f_conv_var_name }} => fortran_type_f%{{ var.name }}
  {% endfor %}
  {% for var in type.variables -%}
  {% if var.is_f_type("Type") -%}
  call {{ var.typedecl.name }}_f_read(lun, {{ var.f_conv_var_name }})
  {% else -%}
  read(UNIT=lun) fortran_type_f%{{ var.name }}
  {% endif -%}
  {% endfor %}
end subroutine {{ type.f_type_name }}_f_read
{% endfor %}

end module {{ module_name }}_io
