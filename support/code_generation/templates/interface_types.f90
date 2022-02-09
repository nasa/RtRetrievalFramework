module {{ module_name }}

use iso_c_binding
{%- for mod_name in addl_used_modules %}
use {{ mod_name }}
{%- endfor %}

! This module was auto-generated 

implicit none

{% if pars_wrapper %}
! Links to module: "{{ pars_wrapper.f_type_name }}" in file: "{{ pars_wrapper.source_basename }}"
type, bind(c) :: {{ pars_wrapper.f_bind_name }}
  {% for var in pars_wrapper.variables -%}
  {%- if var.name == "lidort_version_number" -%}
  {{ var.f_store_var_decls(with_prefix=False, srch_type="argument_char")|f_decl_lines|indent(2) }}({{ var.typedecl.get_length() }})
  {%- else -%}
  {{ var.f_store_var_decls(with_prefix=False)|f_decl_lines|indent(2) }}
  {%- endif %}
  {% endfor %}
end type {{ pars_wrapper.f_bind_name }}
{% endif %}

{% for type in type_classes %}
! Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
type, bind(c) :: {{ type.f_bind_name }}
  {% for var in type.variables -%}
  {% if not var.is_f_type("Character") -%}
  {{ var.f_store_var_decls(with_prefix=False, c_pointer=True)|f_decl_lines|indent(2) }} ! {{ var.f_var_comment() }}
  {%- endif %}
  {{ var.f_supporting_var_decls()|f_decl_lines|indent(2) }}

  {% endfor %}
end type {{ type.f_bind_name }}
{% endfor %}

contains

{% if pars_wrapper %}
subroutine set_lidort_pars(pars_struct) bind(C)
  type({{ pars_wrapper.f_type_name }}_c) :: pars_struct
  integer :: len_idx

  {% for var in pars_wrapper.variables -%}
  {%- if var.name == "lidort_version_number" -%}
  do len_idx = 1, {{var.typedecl.get_length()}}
    pars_struct%{{ var.name }}(len_idx:len_idx) = {{ var.name|upper }}(len_idx:len_idx)
  end do
  {%- else -%}
  pars_struct%{{ var.name }} = {{ var.name|upper }}
  {%- endif %}
  {% endfor %}
end subroutine set_lidort_pars
{% endif %}

{% for type in type_classes %}
! Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
! Allocs and initializes type
subroutine {{ type.f_type_name }}_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  type({{ type.f_type_name }}_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type({{ type.f_type_name }}), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call {{ type.f_type_name }}_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine {{ type.f_type_name }}_c_alloc_init

! Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
! Initializes only with no allocation
subroutine {{ type.f_type_name }}_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use {{ type.f_parent_name }}
  {%- for use_name in type.use_modules %}
  use {{ use_name }}{% endfor %}

  type({{ type.f_type_name }}_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type({{ type.f_type_name }}), pointer :: fortran_type_f

  {% for var in type.variables if not var.is_f_type("Character") -%}
  {% if var.is_f_type("Type") -%}
  type({{ var.typedecl.name }}), pointer :: {{ var.name }}_ptr
  {% else -%}
  {{ var.f_declaration(no_dim_extent=True, f_pointer=True) }} :: {{ var.name+"_ptr" }}
  {% endif -%}
  {% endfor %}

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  {% for var in type.variables -%}
  {% if var.is_f_type("Character") -%}
  fortran_type_f%{{ var.name }} = {{ var.f_zero_value() }}
  transfer_struct_c%{{ var.f_store_var_name }}_f_len = len(fortran_type_f%{{ var.name }})
  {%- else -%}
  {% if not var.is_f_type("Type") -%}{#- Nested type within a type #}
  fortran_type_f%{{ var.name }} = {{ var.f_zero_value() }}
  {% endif %}{#- endif for else on if var.is_f_type("Type") -#}
  {{ var.name }}_ptr => fortran_type_f%{{ var.name }}
  transfer_struct_c%{{ var.f_store_var_name }} = c_loc({{ var.name }}_ptr{{ var.f_lower_bounds() }})
  inquire(iolength=transfer_struct_c%{{ var.name }}_f_byte_size) fortran_type_f%{{ var.name }}{{ var.f_lower_bounds() }}
#ifdef ifort
  transfer_struct_c%{{ var.name }}_f_byte_size = transfer_struct_c%{{ var.name }}_f_byte_size * 4
#endif
  {% endif %}
  {% if var.is_array() %}{% for dim_idx in range(var.get_array_spec()|count) -%}
  transfer_struct_c%{{ var.f_store_var_name }}_f_shapes({{ dim_idx + 1 }}) = size(fortran_type_f%{{ var.name }}, {{ dim_idx + 1}})
  {% endfor %}{% endif %}
  {% endfor %}
end subroutine {{ type.f_type_name }}_c_init_only

subroutine {{ type.f_type_name }}_c_destroy(fortran_type_c) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  type(c_ptr), intent(inout) :: fortran_type_c

  type({{ type.f_type_name }}), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine {{ type.f_type_name }}_c_destroy

subroutine {{ type.f_type_name }}_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type({{ type.f_type_name }}), pointer :: fortran_type_f_from
  type({{ type.f_type_name }}), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  {% for var in type.variables -%}
  fortran_type_f_to%{{ var.name }} = fortran_type_f_from%{{ var.name }}
  {% endfor %}

end subroutine {{ type.f_type_name }}_c_copy
{% endfor %}

{% for type in type_classes %}{% for var in type.variables if var.is_f_type("Character") -%}
! Wrapper for character variable "{{ var.name }}" of type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
subroutine {{ var.f_wrapper_routine_name() }}(fortran_type_c, {{ var.f_arg_var_decls()|f_arg_name_list|indent(6) }}) bind(C)
  use {{ type.f_parent_name }}, only : {{ type.f_type_name }}

  type(c_ptr), intent(inout) :: fortran_type_c
  {{ var.f_arg_var_decls()|f_decl_lines|indent(2) }}

  type({{ type.f_type_name }}), pointer :: fortran_type_f
  {{ var.f_local_var_decls(conv_var=False)|f_decl_lines|indent(2) }}

  call c_f_pointer(fortran_type_c, fortran_type_f)

  {{ var.f_copy_to_c_code(source_name="fortran_type_f%"+var.f_store_var_name, bounds_name=var.f_store_var_name)|join("\n")|indent(2) }}

end subroutine {{ var.f_wrapper_routine_name() }}
{{- type.addl_f_routines|indent(2) }}
{% endfor %}{% endfor %}

end module {{ module_name }}
