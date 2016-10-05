{% from "constructor_def.macro" import constructor -%}
{% from "include_array.macro" import includes -%}
#ifndef {{ file_name|upper }}_H
#define {{ file_name|upper }}_H

#include <iostream>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "fp_exception.h"
{{ includes(addl_includes) }}

/* This file was auto-generated */

#define FORTRAN_TRUE_INT 1

#define BYTE_SIZE_ERROR_CHECK(var_name, c_size, f_size) \
  if(c_size != f_size) { \
    std::stringstream err_msg; \
    err_msg << "Size of C variable: " << c_size \
            << " for " << var_name \
            << " does not match size of Fortran variable: " << f_size; \
    throw Exception(err_msg.str()); \
  }

namespace FullPhysics {

{% if pars_wrapper -%}
extern "C" {
  void set_lidort_pars(struct {{ pars_wrapper.class_name }} *lidort_pars_struct_c);
}

/* This struct cannot inherit any properties from another struct/class which is a requirement of Fortran interoperability */

struct {{ pars_wrapper.class_name }} {

  {% for var in pars_wrapper.variables -%}
  {{ var.c_store_var_decls(with_prefix=False, const=True)|c_decl_lines|indent(2) }}
  {% endfor %}
  static {{ pars_wrapper.class_name }}& instance() {
    static {{ pars_wrapper.class_name }} obj;
    return obj;
  }

  {%set out_op = pars_wrapper.c_output_operator() %}
  {{ out_op.signature }} {
    {{ out_op.contents|indent(4) }}
  }

private:
  {{ constructor(pars_wrapper.class_name, inits=pars_wrapper.class_inits) }} { 
    set_lidort_pars(this);
  }
};
{% endif %}

class {{ type_base_class_name }} : public Printable<{{ type_base_class_name }}> {
public:
  {{ type_base_class_name }}() : fortran_type_c(0), owns_pointer(true) {}
  {{ type_base_class_name }}(void* allocated_f_type_c) : fortran_type_c(allocated_f_type_c), owns_pointer(false) {}
  void* fortran_type_ptr() { return fortran_type_c; }

  virtual void print(std::ostream &output_stream) const = 0;

protected:
  void *fortran_type_c;
  bool owns_pointer;
};

{% for type in type_classes -%}
// Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
extern "C" {
  void {{ type.f_type_name }}_c_alloc_init(struct {{ type.f_type_name }} *transfer_struct_c, void **fortran_type_c);
  void {{ type.f_type_name }}_c_init_only(struct {{ type.f_type_name }} *transfer_struct_c, void **fortran_type_c);
  void {{ type.f_type_name }}_c_copy(void **fortran_type_c_from, void **fortran_type_c_to);
  void {{ type.f_type_name }}_c_destroy(void **fortran_type_c);
  {% for var in type.variables if var.is_f_type("Character") -%}
  void {{ var.f_wrapper_routine_name() }}(void **fortran_type_c, {{ var.c_extern_arg_decls()|c_arg_decl_list }});
  {% endfor %}
}

struct {{ type.f_type_name }} {
  {% for var in type.variables -%}
  {% if var.is_f_type("Type") -%}
  void* {{ var.c_store_var_name }};
  {%- elif not var.is_f_type("Character") -%}
  {{ var.c_store_var_decls(pointer=True, with_prefix=False)|c_decl_lines }}
  {%- endif %}
  {{ var.c_supporting_vars()|join("\n")|indent(2) }}

  {% endfor %}
};

// Links to type: "{{ type.f_type_name }}" from module: "{{ type.f_parent_name }}" in file: "{{ type.source_basename }}"
class {{ type.class_name }} : public {{ type_base_class_name }} {
public:
  // Allocating constructor
  {{ type.class_name }}() : {{ type_base_class_name }}() {
    {{ type.f_type_name }}_c_alloc_init(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Linked to other data structure
  {{ type.class_name }}(void* allocated_f_type_c) :  {{ type_base_class_name }}(allocated_f_type_c) {
    {{ type.f_type_name }}_c_init_only(&transfer_struct_c, &fortran_type_c);
    link_blitz_arrays();
    link_nested_types();
  }

  // Deallocate
  ~{{ type.class_name }}() {
    if (owns_pointer)
      {{ type.f_type_name }}_c_destroy(&fortran_type_c);
  }

  {% for var in type.variables -%}
  {% if var.is_f_type("Type") -%}
  {{ var.c_get_accessor_code(return_const=False, method_const=False)|indent(2) }}
  {{ var.c_get_accessor_code()|indent(2) }}
  {% else -%}
  {{ var.c_get_accessor_code()|indent(2) }}
  {% endif -%}
  {% if not var.is_f_type("Character") -%}
  {{ var.c_set_accessor_code()|indent(2) }}
  {% endif %}
  {% endfor %}

  {%set out_op = type.c_output_operator(use_accessor=True, as_printable=True) %}
  {{ out_op.signature }} {
    {{ out_op.contents|indent(4) }}
  }

  void check_byte_sizes() {
    {% for var in type.variables if not var.is_f_type("Character") and not var.is_f_type("Type") -%}
    BYTE_SIZE_ERROR_CHECK("{{ var.c_store_var_name }}",sizeof(*transfer_struct_c.{{ var.c_store_var_name }}),transfer_struct_c.{{ var.c_store_var_name }}_f_byte_size);
    {% endfor %}
  }
  {{- type.addl_c_public_methods|indent(2) }}

private:
  void link_blitz_arrays() {
    {% for var in type.variables if not var.is_f_type("Character") and var.is_array() -%}
    {{ var.c_store_var_name }}.reference({{ var.blitz_type(use_bool=False) }}(transfer_struct_c.{{ var.c_store_var_name }},
      blitz::shape({%- for dim_idx in range(var.get_array_spec()|count) %}{% if not loop.first %}{{ " "*19 }}{% endif %}transfer_struct_c.{{ var.c_store_var_name }}_f_shapes[{{dim_idx}}]{% if not loop.last %}{{ ",\n" }}{% endif %}{% endfor %}),
      blitz::neverDeleteData, blitz::ColumnMajorArray<{{ var.get_array_spec()|count }}>()));
    {% endfor %}
  }

  void link_nested_types() {
    {% for var in type.variables if var.is_f_type("Type") -%}
      {{ var.c_store_var_name }}.reset( new {{ var.class_name() }}(transfer_struct_c.{{ var.c_store_var_name }}) );
    {% endfor %}
  }

  struct {{ type.f_type_name }} transfer_struct_c;

  {% for var in type.variables if not var.is_f_type("Character") and var.is_array() or var.is_f_type("Type") -%}
  {{ var.c_store_var_decls(blitz_array=True, use_bool=False)|c_decl_lines }}
  {% endfor %}
  {{- type.addl_c_private_methods|indent(2) }}
};

{% endfor %}

}
#endif
