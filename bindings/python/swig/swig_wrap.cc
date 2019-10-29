#include <Python.h>
#include <iostream>
#include "swig_type_mapper_base.h"

// Python 2 and 3 do strings differently, so we have a simple macro to 
// keep from having lots of ifdefs spread around. See 
// https://wiki.python.org/moin/PortingExtensionModulesToPy3k for details on
// this.
#if PY_MAJOR_VERSION > 2
#define Text_FromUTF8(str) PyUnicode_FromString(str)
#else
#define Text_FromUTF8(str) PyString_FromString(str)
#endif

// Python 2 and 3 have different name for their swig init functions
#if PY_MAJOR_VERSION > 2
#define SWIG_INIT_FUNC(S) PyInit__ ## S
#define SWIG_INIT_TYPE PyObject *
#define SWIG_INIT_MODULE init_extension_module3
#else
#define SWIG_INIT_FUNC(S) init_ ## S
#define SWIG_INIT_TYPE void
#define SWIG_INIT_MODULE init_extension_module2
#endif
using namespace FullPhysics;
// Map used between type_index and object to map this to python.
std::map<type_index, boost::shared_ptr<SwigTypeMapperBase> > 
  FullPhysics::swig_type_map;

extern "C" {
#if PY_MAJOR_VERSION > 2
  PyObject * PyInit__swig_wrap(void);
#else
  void init_swig_wrap(void);
#endif
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fp_exception)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(turn_on_fe_exception)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(generic_object)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(unit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(double_with_unit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(constant)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(default_constant)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(array_with_unit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(auto_derivative)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(auto_derivative_with_unit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(array_ad)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(array_ad_with_unit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectral_domain)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(logger)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fp_time)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(heritage_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(linear_algebra)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(rayleigh_greek_moment)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(hdf_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(hdf_sounding_id)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(acos_sounding_id)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_sounding_id)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(uq_sounding_id)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(acos_met_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_met_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(uq_ecmwf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_sim_met_ecmwf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(observer)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fts_run_log)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectral_range)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectral_bound)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectrum)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(named_spectrum)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(state_vector)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(sub_state_vector_array)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(sub_state_vector_proxy)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(hdf_file_generating)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(closest_point)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(polynomial_eval)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(linear_interpolate)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_absorption_spectrum)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_continuum_spectrum)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_doppler_shift)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectrum_effect)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(instrument)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_function)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(instrument_correction)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectrum_effect_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(dispersion)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(initial_guess)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(composite_initial_guess)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(perturbation)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_function)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(convergence_check)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(noise_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_hdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectrum_sampling)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(forward_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(forward_model_spectral_grid)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectral_window)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(stokes_coefficient)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(stokes_coefficient_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(gas_absorption)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature_offset)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(altitude)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_scaled)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ground)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(rt_atmosphere)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiative_transfer)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiative_transfer_retrievable)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiative_transfer_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiative_transfer_fixed_stokes_coefficient)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiative_transfer_single_wn)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spurr_rt)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(hres_wrapper)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spurr_driver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_property)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_property_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_extinction)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_extinction_imp_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_func)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_func_diff)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_problem)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(problem_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_func_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_func_diff_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_problem_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(model_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(iterative_solver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(iterative_solver_der)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_solver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_minimizer)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(model_measure)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(max_likelihood)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(max_a_posteriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure_holder)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(meteorology)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_property_hdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_met_prior)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(co2_profile_prior)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_property_rh_hdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(hdf_constant)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_continuum_table)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_absorption_table)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_doppler_shift_polynomial)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_doppler_shift_l1b)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_absorption_and_continuum)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(bad_sample_noise_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(composite_perturbation)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(dispersion_polynomial)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_table)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_gaussian)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_convolution)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_instrument)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(zero_offset_waveform)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(empirical_orthogonal_function)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiance_scaling)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiance_scaling_sv_fit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiance_scaling_linear_fit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_fts)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ils_fts)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_heritage)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_acos)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_average)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_scale_radiance)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_oco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_uq)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_cache)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_noise_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(gosat_noise_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(precomputed_noise_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absco_hdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(altitude_hydrostatic)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_absco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ground_lambertian)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ground_coxmunk)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ground_coxmunk_plus_lambertian)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(ground_brdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(initial_guess_value)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_optical)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(merra_aerosol)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(rayleigh)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(atmosphere_oco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(relative_humidity)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(connor_solver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(log_timing)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(chisq_convergence)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(connor_convergence)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solver_iteration_log)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(forward_model_cost_function)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(error_analysis)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_forward_model)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(output_hdf)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(uniform_spectrum_sampling)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectrum_sampling_fixed_spacing)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fp_logger)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lidort_interface_types)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lidort_interface_masters)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lidort_driver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lidort_rt)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(twostream_interface)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(twostream_driver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(twostream_rt)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(l_rad_driver)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(l_rad_rt)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lsi_rt)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fm_nlls_problem)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_extinction_linear)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_extinction_log)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(aerosol_shape_gaussian)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(l2_fp_configuration)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(l2_fp_configuration_lua)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nonuniform_spectrum_sampling)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(spectral_window_range)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature_met)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature_level_offset)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_met)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_level)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_log_level)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_level_scaled)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure_sigma)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(stokes_coefficient_constant)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(stokes_coefficient_fraction)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(tccon_apriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(reference_vmr_apriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(gas_vmr_apriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_sim_apriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(dispersion_fit)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fluorescence_effect)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_absorption_gfit_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_solver_gsl)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_solver_gsl_lmsder)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_solver_gsl_lmder)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(cost_minimizer_gsl)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(model_measure_oco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(max_likelihood_oco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(max_a_posteriori_oco)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(connor_solver_map)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_max_likelihood)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_max_a_posteriori)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(nlls_problem_scaled)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure_level_input)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(pressure_fixed_level)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(temperature_fixed_level)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_fixed_level)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(absorber_vmr_fixed_level_scaled)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_continuum_polynomial)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(solar_absorption_oco_file)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(register_output_base)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(dispersion_polynomial_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(stokes_coefficient_fraction_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(zero_offset_waveform_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(empirical_orthogonal_function_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(level_1b_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(state_vector_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(forward_model_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(oco_forward_model_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(connor_convergence_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(connor_solver_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(error_analysis_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(dispersion_fit_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(max_a_posteriori_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(fluorescence_effect_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(radiance_scaling_output)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(lua_state)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(luabind_object)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(swig_std)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(swig_array)(void);
  SWIG_INIT_TYPE SWIG_INIT_FUNC(swig_rational)(void);
}

// Used throughout SWIG wrapper, define here because it is convenient.
std::string parse_python_exception() {
  PyObject *type = NULL, *value = NULL, *tb = NULL;
  PyErr_Fetch(&type, &value, &tb);
  PyObject* mod = PyImport_ImportModule("traceback");
  PyObject* err_str_list = NULL;
  if(tb)
    err_str_list = PyObject_CallMethodObjArgs(mod,
	      Text_FromUTF8("format_exception"), type, value, tb, NULL);
  std::string ret = "Python error that I can't parse";
  if(err_str_list) {
    PyObject* err_str = 
      PyObject_CallMethodObjArgs(Text_FromUTF8(""),
				 Text_FromUTF8("join"), 
				 err_str_list, NULL);
    if(err_str) {
        PyObject * temp_bytes = PyUnicode_AsEncodedString(err_str, "ASCII", 
	"strict");
        ret = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
        Py_DECREF(temp_bytes);
    }
    Py_XDECREF(err_str);
  } else if(value) {
    PyObject * temp_bytes = PyUnicode_AsEncodedString(value, "ASCII", 
	"strict");
    ret = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
    Py_DECREF(temp_bytes);
  }
  Py_XDECREF(mod);
  Py_XDECREF(err_str_list);
  Py_XDECREF(type);
  Py_XDECREF(value);
  Py_XDECREF(tb);
  return ret;
}

#if PY_MAJOR_VERSION > 2
// Version for python 3
static void init_extension_module3(PyObject* package, const char *modulename,
				  PyObject * (*initfunction)(void)) {
  PyObject *module = initfunction();
  PyObject *module_dic = PyImport_GetModuleDict();
  PyDict_SetItem(module_dic, Text_FromUTF8(modulename), module);
  if(PyModule_AddObject(package, (char *)modulename, module)) {
    std::cerr << "Initialisation in PyImport_AddObject failed for module "
	      << modulename << "\n";
    return;
  }
  Py_INCREF(module);
}
#else 
// Version for python 2
static void init_extension_module2(PyObject* package, const char *modulename,
				  void (*initfunction)(void)) {
  PyObject *module = PyImport_AddModule((char *)modulename);
  if(!module) {
    std::cerr << "Initialisation in PyImport_AddModule failed for module "
	      << modulename << "\n";
    return;
  }
  if(PyModule_AddObject(package, (char *)modulename, module)) {
    std::cerr << "Initialisation in PyImport_AddObject failed for module "
	      << modulename << "\n";
    return;
  }
  Py_INCREF(module);
  initfunction();
}
#endif


// This next blob of code comes from 
// https://wiki.python.org/moin/PortingExtensionModulesToPy3k

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef swig_wrap_methods[] = {
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int swig_wrap_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int swig_wrap_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_swig_wrap",
        NULL,
        sizeof(struct module_state),
        swig_wrap_methods,
        NULL,
        swig_wrap_traverse,
        swig_wrap_clear,
        NULL
};

#define INITERROR return NULL

PyObject *
PyInit__swig_wrap(void)

#else
#define INITERROR return

void
init_swig_wrap(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("_swig_wrap", swig_wrap_methods);
#endif

    if (module == NULL) {
        std::cerr << "Initialization failed\n";
        INITERROR;
    }
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("swig_wrap.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

  PyObject *package = PyImport_AddModule((char *)"new_geocal");
  if(!package) {
      std::cerr << "Initialization failed\n";
      INITERROR;
  }
  
  SWIG_INIT_MODULE(package, "_fp_exception", SWIG_INIT_FUNC(fp_exception));
  SWIG_INIT_MODULE(package, "_turn_on_fe_exception", SWIG_INIT_FUNC(turn_on_fe_exception));
  SWIG_INIT_MODULE(package, "_generic_object", SWIG_INIT_FUNC(generic_object));
  SWIG_INIT_MODULE(package, "_unit", SWIG_INIT_FUNC(unit));
  SWIG_INIT_MODULE(package, "_double_with_unit", SWIG_INIT_FUNC(double_with_unit));
  SWIG_INIT_MODULE(package, "_constant", SWIG_INIT_FUNC(constant));
  SWIG_INIT_MODULE(package, "_default_constant", SWIG_INIT_FUNC(default_constant));
  SWIG_INIT_MODULE(package, "_array_with_unit", SWIG_INIT_FUNC(array_with_unit));
  SWIG_INIT_MODULE(package, "_auto_derivative", SWIG_INIT_FUNC(auto_derivative));
  SWIG_INIT_MODULE(package, "_auto_derivative_with_unit", SWIG_INIT_FUNC(auto_derivative_with_unit));
  SWIG_INIT_MODULE(package, "_array_ad", SWIG_INIT_FUNC(array_ad));
  SWIG_INIT_MODULE(package, "_array_ad_with_unit", SWIG_INIT_FUNC(array_ad_with_unit));
  SWIG_INIT_MODULE(package, "_spectral_domain", SWIG_INIT_FUNC(spectral_domain));
  SWIG_INIT_MODULE(package, "_logger", SWIG_INIT_FUNC(logger));
  SWIG_INIT_MODULE(package, "_fp_time", SWIG_INIT_FUNC(fp_time));
  SWIG_INIT_MODULE(package, "_heritage_file", SWIG_INIT_FUNC(heritage_file));
  SWIG_INIT_MODULE(package, "_linear_algebra", SWIG_INIT_FUNC(linear_algebra));
  SWIG_INIT_MODULE(package, "_rayleigh_greek_moment", SWIG_INIT_FUNC(rayleigh_greek_moment));
  SWIG_INIT_MODULE(package, "_hdf_file", SWIG_INIT_FUNC(hdf_file));
  SWIG_INIT_MODULE(package, "_hdf_sounding_id", SWIG_INIT_FUNC(hdf_sounding_id));
  SWIG_INIT_MODULE(package, "_acos_sounding_id", SWIG_INIT_FUNC(acos_sounding_id));
  SWIG_INIT_MODULE(package, "_oco_sounding_id", SWIG_INIT_FUNC(oco_sounding_id));
  SWIG_INIT_MODULE(package, "_uq_sounding_id", SWIG_INIT_FUNC(uq_sounding_id));
  SWIG_INIT_MODULE(package, "_acos_met_file", SWIG_INIT_FUNC(acos_met_file));
  SWIG_INIT_MODULE(package, "_oco_met_file", SWIG_INIT_FUNC(oco_met_file));
  SWIG_INIT_MODULE(package, "_uq_ecmwf", SWIG_INIT_FUNC(uq_ecmwf));
  SWIG_INIT_MODULE(package, "_oco_sim_met_ecmwf", SWIG_INIT_FUNC(oco_sim_met_ecmwf));
  SWIG_INIT_MODULE(package, "_observer", SWIG_INIT_FUNC(observer));
  SWIG_INIT_MODULE(package, "_fts_run_log", SWIG_INIT_FUNC(fts_run_log));
  SWIG_INIT_MODULE(package, "_spectral_range", SWIG_INIT_FUNC(spectral_range));
  SWIG_INIT_MODULE(package, "_spectral_bound", SWIG_INIT_FUNC(spectral_bound));
  SWIG_INIT_MODULE(package, "_spectrum", SWIG_INIT_FUNC(spectrum));
  SWIG_INIT_MODULE(package, "_named_spectrum", SWIG_INIT_FUNC(named_spectrum));
  SWIG_INIT_MODULE(package, "_state_vector", SWIG_INIT_FUNC(state_vector));
  SWIG_INIT_MODULE(package, "_sub_state_vector_array", SWIG_INIT_FUNC(sub_state_vector_array));
  SWIG_INIT_MODULE(package, "_sub_state_vector_proxy", SWIG_INIT_FUNC(sub_state_vector_proxy));
  SWIG_INIT_MODULE(package, "_hdf_file_generating", SWIG_INIT_FUNC(hdf_file_generating));
  SWIG_INIT_MODULE(package, "_closest_point", SWIG_INIT_FUNC(closest_point));
  SWIG_INIT_MODULE(package, "_polynomial_eval", SWIG_INIT_FUNC(polynomial_eval));
  SWIG_INIT_MODULE(package, "_linear_interpolate", SWIG_INIT_FUNC(linear_interpolate));
  SWIG_INIT_MODULE(package, "_solar_absorption_spectrum", SWIG_INIT_FUNC(solar_absorption_spectrum));
  SWIG_INIT_MODULE(package, "_solar_continuum_spectrum", SWIG_INIT_FUNC(solar_continuum_spectrum));
  SWIG_INIT_MODULE(package, "_solar_doppler_shift", SWIG_INIT_FUNC(solar_doppler_shift));
  SWIG_INIT_MODULE(package, "_spectrum_effect", SWIG_INIT_FUNC(spectrum_effect));
  SWIG_INIT_MODULE(package, "_solar_model", SWIG_INIT_FUNC(solar_model));
  SWIG_INIT_MODULE(package, "_instrument", SWIG_INIT_FUNC(instrument));
  SWIG_INIT_MODULE(package, "_aerosol", SWIG_INIT_FUNC(aerosol));
  SWIG_INIT_MODULE(package, "_ils", SWIG_INIT_FUNC(ils));
  SWIG_INIT_MODULE(package, "_ils_function", SWIG_INIT_FUNC(ils_function));
  SWIG_INIT_MODULE(package, "_instrument_correction", SWIG_INIT_FUNC(instrument_correction));
  SWIG_INIT_MODULE(package, "_spectrum_effect_imp_base", SWIG_INIT_FUNC(spectrum_effect_imp_base));
  SWIG_INIT_MODULE(package, "_dispersion", SWIG_INIT_FUNC(dispersion));
  SWIG_INIT_MODULE(package, "_initial_guess", SWIG_INIT_FUNC(initial_guess));
  SWIG_INIT_MODULE(package, "_composite_initial_guess", SWIG_INIT_FUNC(composite_initial_guess));
  SWIG_INIT_MODULE(package, "_perturbation", SWIG_INIT_FUNC(perturbation));
  SWIG_INIT_MODULE(package, "_cost_function", SWIG_INIT_FUNC(cost_function));
  SWIG_INIT_MODULE(package, "_convergence_check", SWIG_INIT_FUNC(convergence_check));
  SWIG_INIT_MODULE(package, "_level_1b", SWIG_INIT_FUNC(level_1b));
  SWIG_INIT_MODULE(package, "_noise_model", SWIG_INIT_FUNC(noise_model));
  SWIG_INIT_MODULE(package, "_level_1b_hdf", SWIG_INIT_FUNC(level_1b_hdf));
  SWIG_INIT_MODULE(package, "_spectrum_sampling", SWIG_INIT_FUNC(spectrum_sampling));
  SWIG_INIT_MODULE(package, "_forward_model", SWIG_INIT_FUNC(forward_model));
  SWIG_INIT_MODULE(package, "_forward_model_spectral_grid", SWIG_INIT_FUNC(forward_model_spectral_grid));
  SWIG_INIT_MODULE(package, "_output", SWIG_INIT_FUNC(output));
  SWIG_INIT_MODULE(package, "_spectral_window", SWIG_INIT_FUNC(spectral_window));
  SWIG_INIT_MODULE(package, "_pressure", SWIG_INIT_FUNC(pressure));
  SWIG_INIT_MODULE(package, "_stokes_coefficient", SWIG_INIT_FUNC(stokes_coefficient));
  SWIG_INIT_MODULE(package, "_pressure_imp_base", SWIG_INIT_FUNC(pressure_imp_base));
  SWIG_INIT_MODULE(package, "_stokes_coefficient_imp_base", SWIG_INIT_FUNC(stokes_coefficient_imp_base));
  SWIG_INIT_MODULE(package, "_gas_absorption", SWIG_INIT_FUNC(gas_absorption));
  SWIG_INIT_MODULE(package, "_temperature", SWIG_INIT_FUNC(temperature));
  SWIG_INIT_MODULE(package, "_temperature_imp_base", SWIG_INIT_FUNC(temperature_imp_base));
  SWIG_INIT_MODULE(package, "_temperature_offset", SWIG_INIT_FUNC(temperature_offset));
  SWIG_INIT_MODULE(package, "_altitude", SWIG_INIT_FUNC(altitude));
  SWIG_INIT_MODULE(package, "_absorber_vmr", SWIG_INIT_FUNC(absorber_vmr));
  SWIG_INIT_MODULE(package, "_absorber_vmr_imp_base", SWIG_INIT_FUNC(absorber_vmr_imp_base));
  SWIG_INIT_MODULE(package, "_absorber_vmr_scaled", SWIG_INIT_FUNC(absorber_vmr_scaled));
  SWIG_INIT_MODULE(package, "_absorber", SWIG_INIT_FUNC(absorber));
  SWIG_INIT_MODULE(package, "_ground", SWIG_INIT_FUNC(ground));
  SWIG_INIT_MODULE(package, "_rt_atmosphere", SWIG_INIT_FUNC(rt_atmosphere));
  SWIG_INIT_MODULE(package, "_radiative_transfer", SWIG_INIT_FUNC(radiative_transfer));
  SWIG_INIT_MODULE(package, "_radiative_transfer_retrievable", SWIG_INIT_FUNC(radiative_transfer_retrievable));
  SWIG_INIT_MODULE(package, "_radiative_transfer_imp_base", SWIG_INIT_FUNC(radiative_transfer_imp_base));
  SWIG_INIT_MODULE(package, "_radiative_transfer_fixed_stokes_coefficient", SWIG_INIT_FUNC(radiative_transfer_fixed_stokes_coefficient));
  SWIG_INIT_MODULE(package, "_radiative_transfer_single_wn", SWIG_INIT_FUNC(radiative_transfer_single_wn));
  SWIG_INIT_MODULE(package, "_spurr_rt", SWIG_INIT_FUNC(spurr_rt));
  SWIG_INIT_MODULE(package, "_hres_wrapper", SWIG_INIT_FUNC(hres_wrapper));
  SWIG_INIT_MODULE(package, "_spurr_driver", SWIG_INIT_FUNC(spurr_driver));
  SWIG_INIT_MODULE(package, "_aerosol_property", SWIG_INIT_FUNC(aerosol_property));
  SWIG_INIT_MODULE(package, "_aerosol_property_imp_base", SWIG_INIT_FUNC(aerosol_property_imp_base));
  SWIG_INIT_MODULE(package, "_aerosol_extinction", SWIG_INIT_FUNC(aerosol_extinction));
  SWIG_INIT_MODULE(package, "_aerosol_extinction_imp_base", SWIG_INIT_FUNC(aerosol_extinction_imp_base));
  SWIG_INIT_MODULE(package, "_cost_func", SWIG_INIT_FUNC(cost_func));
  SWIG_INIT_MODULE(package, "_cost_func_diff", SWIG_INIT_FUNC(cost_func_diff));
  SWIG_INIT_MODULE(package, "_nlls_problem", SWIG_INIT_FUNC(nlls_problem));
  SWIG_INIT_MODULE(package, "_problem_state", SWIG_INIT_FUNC(problem_state));
  SWIG_INIT_MODULE(package, "_cost_func_state", SWIG_INIT_FUNC(cost_func_state));
  SWIG_INIT_MODULE(package, "_cost_func_diff_state", SWIG_INIT_FUNC(cost_func_diff_state));
  SWIG_INIT_MODULE(package, "_nlls_problem_state", SWIG_INIT_FUNC(nlls_problem_state));
  SWIG_INIT_MODULE(package, "_model_state", SWIG_INIT_FUNC(model_state));
  SWIG_INIT_MODULE(package, "_iterative_solver", SWIG_INIT_FUNC(iterative_solver));
  SWIG_INIT_MODULE(package, "_iterative_solver_der", SWIG_INIT_FUNC(iterative_solver_der));
  SWIG_INIT_MODULE(package, "_nlls_solver", SWIG_INIT_FUNC(nlls_solver));
  SWIG_INIT_MODULE(package, "_cost_minimizer", SWIG_INIT_FUNC(cost_minimizer));
  SWIG_INIT_MODULE(package, "_model_measure", SWIG_INIT_FUNC(model_measure));
  SWIG_INIT_MODULE(package, "_max_likelihood", SWIG_INIT_FUNC(max_likelihood));
  SWIG_INIT_MODULE(package, "_max_a_posteriori", SWIG_INIT_FUNC(max_a_posteriori));
  SWIG_INIT_MODULE(package, "_pressure_holder", SWIG_INIT_FUNC(pressure_holder));
  SWIG_INIT_MODULE(package, "_meteorology", SWIG_INIT_FUNC(meteorology));
  SWIG_INIT_MODULE(package, "_aerosol_property_hdf", SWIG_INIT_FUNC(aerosol_property_hdf));
  SWIG_INIT_MODULE(package, "_aerosol_met_prior", SWIG_INIT_FUNC(aerosol_met_prior));
  SWIG_INIT_MODULE(package, "_co2_profile_prior", SWIG_INIT_FUNC(co2_profile_prior));
  SWIG_INIT_MODULE(package, "_aerosol_property_rh_hdf", SWIG_INIT_FUNC(aerosol_property_rh_hdf));
  SWIG_INIT_MODULE(package, "_hdf_constant", SWIG_INIT_FUNC(hdf_constant));
  SWIG_INIT_MODULE(package, "_solar_continuum_table", SWIG_INIT_FUNC(solar_continuum_table));
  SWIG_INIT_MODULE(package, "_solar_absorption_table", SWIG_INIT_FUNC(solar_absorption_table));
  SWIG_INIT_MODULE(package, "_solar_doppler_shift_polynomial", SWIG_INIT_FUNC(solar_doppler_shift_polynomial));
  SWIG_INIT_MODULE(package, "_solar_doppler_shift_l1b", SWIG_INIT_FUNC(solar_doppler_shift_l1b));
  SWIG_INIT_MODULE(package, "_solar_absorption_and_continuum", SWIG_INIT_FUNC(solar_absorption_and_continuum));
  SWIG_INIT_MODULE(package, "_bad_sample_noise_model", SWIG_INIT_FUNC(bad_sample_noise_model));
  SWIG_INIT_MODULE(package, "_composite_perturbation", SWIG_INIT_FUNC(composite_perturbation));
  SWIG_INIT_MODULE(package, "_dispersion_polynomial", SWIG_INIT_FUNC(dispersion_polynomial));
  SWIG_INIT_MODULE(package, "_ils_table", SWIG_INIT_FUNC(ils_table));
  SWIG_INIT_MODULE(package, "_ils_gaussian", SWIG_INIT_FUNC(ils_gaussian));
  SWIG_INIT_MODULE(package, "_ils_convolution", SWIG_INIT_FUNC(ils_convolution));
  SWIG_INIT_MODULE(package, "_ils_instrument", SWIG_INIT_FUNC(ils_instrument));
  SWIG_INIT_MODULE(package, "_zero_offset_waveform", SWIG_INIT_FUNC(zero_offset_waveform));
  SWIG_INIT_MODULE(package, "_empirical_orthogonal_function", SWIG_INIT_FUNC(empirical_orthogonal_function));
  SWIG_INIT_MODULE(package, "_radiance_scaling", SWIG_INIT_FUNC(radiance_scaling));
  SWIG_INIT_MODULE(package, "_radiance_scaling_sv_fit", SWIG_INIT_FUNC(radiance_scaling_sv_fit));
  SWIG_INIT_MODULE(package, "_radiance_scaling_linear_fit", SWIG_INIT_FUNC(radiance_scaling_linear_fit));
  SWIG_INIT_MODULE(package, "_level_1b_fts", SWIG_INIT_FUNC(level_1b_fts));
  SWIG_INIT_MODULE(package, "_ils_fts", SWIG_INIT_FUNC(ils_fts));
  SWIG_INIT_MODULE(package, "_level_1b_heritage", SWIG_INIT_FUNC(level_1b_heritage));
  SWIG_INIT_MODULE(package, "_level_1b_acos", SWIG_INIT_FUNC(level_1b_acos));
  SWIG_INIT_MODULE(package, "_level_1b_average", SWIG_INIT_FUNC(level_1b_average));
  SWIG_INIT_MODULE(package, "_level_1b_scale_radiance", SWIG_INIT_FUNC(level_1b_scale_radiance));
  SWIG_INIT_MODULE(package, "_level_1b_oco", SWIG_INIT_FUNC(level_1b_oco));
  SWIG_INIT_MODULE(package, "_level_1b_uq", SWIG_INIT_FUNC(level_1b_uq));
  SWIG_INIT_MODULE(package, "_level_1b_cache", SWIG_INIT_FUNC(level_1b_cache));
  SWIG_INIT_MODULE(package, "_oco_noise_model", SWIG_INIT_FUNC(oco_noise_model));
  SWIG_INIT_MODULE(package, "_gosat_noise_model", SWIG_INIT_FUNC(gosat_noise_model));
  SWIG_INIT_MODULE(package, "_precomputed_noise_model", SWIG_INIT_FUNC(precomputed_noise_model));
  SWIG_INIT_MODULE(package, "_absco", SWIG_INIT_FUNC(absco));
  SWIG_INIT_MODULE(package, "_absco_hdf", SWIG_INIT_FUNC(absco_hdf));
  SWIG_INIT_MODULE(package, "_altitude_hydrostatic", SWIG_INIT_FUNC(altitude_hydrostatic));
  SWIG_INIT_MODULE(package, "_absorber_absco", SWIG_INIT_FUNC(absorber_absco));
  SWIG_INIT_MODULE(package, "_ground_lambertian", SWIG_INIT_FUNC(ground_lambertian));
  SWIG_INIT_MODULE(package, "_ground_coxmunk", SWIG_INIT_FUNC(ground_coxmunk));
  SWIG_INIT_MODULE(package, "_ground_coxmunk_plus_lambertian", SWIG_INIT_FUNC(ground_coxmunk_plus_lambertian));
  SWIG_INIT_MODULE(package, "_ground_brdf", SWIG_INIT_FUNC(ground_brdf));
  SWIG_INIT_MODULE(package, "_initial_guess_value", SWIG_INIT_FUNC(initial_guess_value));
  SWIG_INIT_MODULE(package, "_aerosol_optical", SWIG_INIT_FUNC(aerosol_optical));
  SWIG_INIT_MODULE(package, "_merra_aerosol", SWIG_INIT_FUNC(merra_aerosol));
  SWIG_INIT_MODULE(package, "_rayleigh", SWIG_INIT_FUNC(rayleigh));
  SWIG_INIT_MODULE(package, "_atmosphere_oco", SWIG_INIT_FUNC(atmosphere_oco));
  SWIG_INIT_MODULE(package, "_relative_humidity", SWIG_INIT_FUNC(relative_humidity));
  SWIG_INIT_MODULE(package, "_connor_solver", SWIG_INIT_FUNC(connor_solver));
  SWIG_INIT_MODULE(package, "_log_timing", SWIG_INIT_FUNC(log_timing));
  SWIG_INIT_MODULE(package, "_chisq_convergence", SWIG_INIT_FUNC(chisq_convergence));
  SWIG_INIT_MODULE(package, "_connor_convergence", SWIG_INIT_FUNC(connor_convergence));
  SWIG_INIT_MODULE(package, "_solver_iteration_log", SWIG_INIT_FUNC(solver_iteration_log));
  SWIG_INIT_MODULE(package, "_forward_model_cost_function", SWIG_INIT_FUNC(forward_model_cost_function));
  SWIG_INIT_MODULE(package, "_error_analysis", SWIG_INIT_FUNC(error_analysis));
  SWIG_INIT_MODULE(package, "_oco_forward_model", SWIG_INIT_FUNC(oco_forward_model));
  SWIG_INIT_MODULE(package, "_output_hdf", SWIG_INIT_FUNC(output_hdf));
  SWIG_INIT_MODULE(package, "_uniform_spectrum_sampling", SWIG_INIT_FUNC(uniform_spectrum_sampling));
  SWIG_INIT_MODULE(package, "_spectrum_sampling_fixed_spacing", SWIG_INIT_FUNC(spectrum_sampling_fixed_spacing));
  SWIG_INIT_MODULE(package, "_fp_logger", SWIG_INIT_FUNC(fp_logger));
  SWIG_INIT_MODULE(package, "_lidort_interface_types", SWIG_INIT_FUNC(lidort_interface_types));
  SWIG_INIT_MODULE(package, "_lidort_interface_masters", SWIG_INIT_FUNC(lidort_interface_masters));
  SWIG_INIT_MODULE(package, "_lidort_driver", SWIG_INIT_FUNC(lidort_driver));
  SWIG_INIT_MODULE(package, "_lidort_rt", SWIG_INIT_FUNC(lidort_rt));
  SWIG_INIT_MODULE(package, "_twostream_interface", SWIG_INIT_FUNC(twostream_interface));
  SWIG_INIT_MODULE(package, "_twostream_driver", SWIG_INIT_FUNC(twostream_driver));
  SWIG_INIT_MODULE(package, "_twostream_rt", SWIG_INIT_FUNC(twostream_rt));
  SWIG_INIT_MODULE(package, "_l_rad_driver", SWIG_INIT_FUNC(l_rad_driver));
  SWIG_INIT_MODULE(package, "_l_rad_rt", SWIG_INIT_FUNC(l_rad_rt));
  SWIG_INIT_MODULE(package, "_lsi_rt", SWIG_INIT_FUNC(lsi_rt));
  SWIG_INIT_MODULE(package, "_fm_nlls_problem", SWIG_INIT_FUNC(fm_nlls_problem));
  SWIG_INIT_MODULE(package, "_aerosol_extinction_linear", SWIG_INIT_FUNC(aerosol_extinction_linear));
  SWIG_INIT_MODULE(package, "_aerosol_extinction_log", SWIG_INIT_FUNC(aerosol_extinction_log));
  SWIG_INIT_MODULE(package, "_aerosol_shape_gaussian", SWIG_INIT_FUNC(aerosol_shape_gaussian));
  SWIG_INIT_MODULE(package, "_l2_fp_configuration", SWIG_INIT_FUNC(l2_fp_configuration));
  SWIG_INIT_MODULE(package, "_l2_fp_configuration_lua", SWIG_INIT_FUNC(l2_fp_configuration_lua));
  SWIG_INIT_MODULE(package, "_nonuniform_spectrum_sampling", SWIG_INIT_FUNC(nonuniform_spectrum_sampling));
  SWIG_INIT_MODULE(package, "_spectral_window_range", SWIG_INIT_FUNC(spectral_window_range));
  SWIG_INIT_MODULE(package, "_temperature_met", SWIG_INIT_FUNC(temperature_met));
  SWIG_INIT_MODULE(package, "_temperature_level_offset", SWIG_INIT_FUNC(temperature_level_offset));
  SWIG_INIT_MODULE(package, "_absorber_vmr_met", SWIG_INIT_FUNC(absorber_vmr_met));
  SWIG_INIT_MODULE(package, "_absorber_vmr_level", SWIG_INIT_FUNC(absorber_vmr_level));
  SWIG_INIT_MODULE(package, "_absorber_vmr_log_level", SWIG_INIT_FUNC(absorber_vmr_log_level));
  SWIG_INIT_MODULE(package, "_absorber_vmr_level_scaled", SWIG_INIT_FUNC(absorber_vmr_level_scaled));
  SWIG_INIT_MODULE(package, "_pressure_sigma", SWIG_INIT_FUNC(pressure_sigma));
  SWIG_INIT_MODULE(package, "_stokes_coefficient_constant", SWIG_INIT_FUNC(stokes_coefficient_constant));
  SWIG_INIT_MODULE(package, "_stokes_coefficient_fraction", SWIG_INIT_FUNC(stokes_coefficient_fraction));
  SWIG_INIT_MODULE(package, "_tccon_apriori", SWIG_INIT_FUNC(tccon_apriori));
  SWIG_INIT_MODULE(package, "_reference_vmr_apriori", SWIG_INIT_FUNC(reference_vmr_apriori));
  SWIG_INIT_MODULE(package, "_gas_vmr_apriori", SWIG_INIT_FUNC(gas_vmr_apriori));
  SWIG_INIT_MODULE(package, "_oco_sim_apriori", SWIG_INIT_FUNC(oco_sim_apriori));
  SWIG_INIT_MODULE(package, "_dispersion_fit", SWIG_INIT_FUNC(dispersion_fit));
  SWIG_INIT_MODULE(package, "_fluorescence_effect", SWIG_INIT_FUNC(fluorescence_effect));
  SWIG_INIT_MODULE(package, "_solar_absorption_gfit_file", SWIG_INIT_FUNC(solar_absorption_gfit_file));
  SWIG_INIT_MODULE(package, "_nlls_solver_gsl", SWIG_INIT_FUNC(nlls_solver_gsl));
  SWIG_INIT_MODULE(package, "_nlls_solver_gsl_lmsder", SWIG_INIT_FUNC(nlls_solver_gsl_lmsder));
  SWIG_INIT_MODULE(package, "_nlls_solver_gsl_lmder", SWIG_INIT_FUNC(nlls_solver_gsl_lmder));
  SWIG_INIT_MODULE(package, "_cost_minimizer_gsl", SWIG_INIT_FUNC(cost_minimizer_gsl));
  SWIG_INIT_MODULE(package, "_model_measure_oco", SWIG_INIT_FUNC(model_measure_oco));
  SWIG_INIT_MODULE(package, "_max_likelihood_oco", SWIG_INIT_FUNC(max_likelihood_oco));
  SWIG_INIT_MODULE(package, "_max_a_posteriori_oco", SWIG_INIT_FUNC(max_a_posteriori_oco));
  SWIG_INIT_MODULE(package, "_connor_solver_map", SWIG_INIT_FUNC(connor_solver_map));
  SWIG_INIT_MODULE(package, "_nlls_max_likelihood", SWIG_INIT_FUNC(nlls_max_likelihood));
  SWIG_INIT_MODULE(package, "_nlls_max_a_posteriori", SWIG_INIT_FUNC(nlls_max_a_posteriori));
  SWIG_INIT_MODULE(package, "_nlls_problem_scaled", SWIG_INIT_FUNC(nlls_problem_scaled));
  SWIG_INIT_MODULE(package, "_pressure_level_input", SWIG_INIT_FUNC(pressure_level_input));
  SWIG_INIT_MODULE(package, "_pressure_fixed_level", SWIG_INIT_FUNC(pressure_fixed_level));
  SWIG_INIT_MODULE(package, "_temperature_fixed_level", SWIG_INIT_FUNC(temperature_fixed_level));
  SWIG_INIT_MODULE(package, "_absorber_vmr_fixed_level", SWIG_INIT_FUNC(absorber_vmr_fixed_level));
  SWIG_INIT_MODULE(package, "_absorber_vmr_fixed_level_scaled", SWIG_INIT_FUNC(absorber_vmr_fixed_level_scaled));
  SWIG_INIT_MODULE(package, "_solar_continuum_polynomial", SWIG_INIT_FUNC(solar_continuum_polynomial));
  SWIG_INIT_MODULE(package, "_solar_absorption_oco_file", SWIG_INIT_FUNC(solar_absorption_oco_file));
  SWIG_INIT_MODULE(package, "_register_output_base", SWIG_INIT_FUNC(register_output_base));
  SWIG_INIT_MODULE(package, "_dispersion_polynomial_output", SWIG_INIT_FUNC(dispersion_polynomial_output));
  SWIG_INIT_MODULE(package, "_stokes_coefficient_fraction_output", SWIG_INIT_FUNC(stokes_coefficient_fraction_output));
  SWIG_INIT_MODULE(package, "_zero_offset_waveform_output", SWIG_INIT_FUNC(zero_offset_waveform_output));
  SWIG_INIT_MODULE(package, "_empirical_orthogonal_function_output", SWIG_INIT_FUNC(empirical_orthogonal_function_output));
  SWIG_INIT_MODULE(package, "_level_1b_output", SWIG_INIT_FUNC(level_1b_output));
  SWIG_INIT_MODULE(package, "_state_vector_output", SWIG_INIT_FUNC(state_vector_output));
  SWIG_INIT_MODULE(package, "_forward_model_output", SWIG_INIT_FUNC(forward_model_output));
  SWIG_INIT_MODULE(package, "_oco_forward_model_output", SWIG_INIT_FUNC(oco_forward_model_output));
  SWIG_INIT_MODULE(package, "_connor_convergence_output", SWIG_INIT_FUNC(connor_convergence_output));
  SWIG_INIT_MODULE(package, "_connor_solver_output", SWIG_INIT_FUNC(connor_solver_output));
  SWIG_INIT_MODULE(package, "_error_analysis_output", SWIG_INIT_FUNC(error_analysis_output));
  SWIG_INIT_MODULE(package, "_dispersion_fit_output", SWIG_INIT_FUNC(dispersion_fit_output));
  SWIG_INIT_MODULE(package, "_max_a_posteriori_output", SWIG_INIT_FUNC(max_a_posteriori_output));
  SWIG_INIT_MODULE(package, "_fluorescence_effect_output", SWIG_INIT_FUNC(fluorescence_effect_output));
  SWIG_INIT_MODULE(package, "_radiance_scaling_output", SWIG_INIT_FUNC(radiance_scaling_output));
  SWIG_INIT_MODULE(package, "_lua_state", SWIG_INIT_FUNC(lua_state));
  SWIG_INIT_MODULE(package, "_luabind_object", SWIG_INIT_FUNC(luabind_object));
  SWIG_INIT_MODULE(package, "_swig_std", SWIG_INIT_FUNC(swig_std));
  SWIG_INIT_MODULE(package, "_swig_array", SWIG_INIT_FUNC(swig_array));
  SWIG_INIT_MODULE(package, "_swig_rational", SWIG_INIT_FUNC(swig_rational));

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


