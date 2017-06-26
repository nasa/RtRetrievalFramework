{% extends "interface_masters.h" %}
{%- block post_wrap %}
{{- super() }}
{%- if wrap.c_routine_name == "run" and master.name.find("brdf") < 0  %}
    Lidort_Pars lid_pars = Lidort_Pars::instance();
    if( lidort_out().status().ts_status_inputcheck() != lid_pars.lidort_success ||
        lidort_out().status().ts_status_calculation() != lid_pars.lidort_success ) {
       std::stringstream err_msg;
       err_msg << "LIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << lidort_out().status();
       throw Exception(err_msg.str());
    }
{%- endif %}
{% endblock %}
