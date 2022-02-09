{% set TEST_TOL = "1e-10" -%}
{% from "include_array.macro" import includes -%}
#include <blitz/array.h>
#include <cmath>
#include "unit_test_support.h"
{{ includes(addl_includes) }}

using namespace FullPhysics;
using namespace blitz;

/* This file was auto-generated */

BOOST_FIXTURE_TEST_SUITE({{ suite_name }}, GlobalFixture)

{% if pars_wrapper -%}
{% set pars_class_name = pars_wrapper.class_name %}
BOOST_AUTO_TEST_CASE({{ pars_class_name|lower }})
{
  // Test obtaining instance
  {{ pars_wrapper.class_name }} tst_obj = {{ pars_wrapper.class_name }}::instance();
  // Test that the value obtained from fortran is as expected
  {% for var in pars_wrapper.variables -%}
  {% if var.is_f_type("Real") -%}
  BOOST_CHECK_CLOSE(tst_obj.{{var.name}}, {{ var.init_string("tst_obj", pars_wrapper.variables) }}, {{ TEST_TOL }});
  {%- else -%}
  BOOST_CHECK_EQUAL(tst_obj.{{var.name}}, {{ var.init_string("tst_obj", pars_wrapper.variables) }});
  {%- endif %}
  {% endfor %}
}
{%- endif %}

{% for type in type_classes -%}
BOOST_AUTO_TEST_CASE({{ type.f_type_name }})
{
  {% if pars_wrapper -%}
  // Used for checking dimensions
  {{ pars_class_name }} lid_pars = {{ pars_class_name }}::instance();
  {%- endif %}

  // Test constructor
  {{ type.class_name }} tst_obj = {{ type.class_name }}();

  // Test memory space sizes
  tst_obj.check_byte_sizes(); // throws exception on mismatch

  // Test variable shapes
  {% for var in type.variables if not var.is_f_type("Character") and var.is_array() -%}
  {% for dim_size in var.expected_dims() -%}
  BOOST_CHECK_EQUAL(tst_obj.{{ var.name }}().extent({{ loop.index0 }}), lid_pars.{{ dim_size }});
  {% endfor -%}
  {%- endfor %}

  // Test initialization
  {% for var in type.variables if not var.is_f_type("Character") and not var.is_f_type("Type") -%}
  {% if var.is_array() -%}
  {{ var.blitz_type(use_bool=True)|trim }} {{ var.name }}_exp(tst_obj.{{ var.name }}().shape());
  {{ var.name }}_exp = {{ var.c_zero_value() }};
  BOOST_CHECK_MATRIX_CLOSE_TOL(tst_obj.{{ var.name }}(), {{ var.name }}_exp, {{ TEST_TOL }});
  {%- else -%}
  {% if var.is_f_type("Real") -%}
  BOOST_CHECK_CLOSE(tst_obj.{{ var.name }}(), {{ var.c_zero_value() }}, {{ TEST_TOL }});
  {%- else -%}
  BOOST_CHECK_EQUAL(tst_obj.{{ var.name }}(), {{ var.c_zero_value() }});
  {%- endif -%}{#- end typeReal #}
  {%- endif -%}{#- not array #}
  {% endfor %}
}

{% endfor %}

BOOST_AUTO_TEST_SUITE_END()
       
