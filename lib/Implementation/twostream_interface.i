// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// Don't bother moving the various function to python attributes. This file is
// *long*, and auto-generated. Not sure that we will use this much directly
// in python.  If this ever becomes an issue we can revisit this, but for
// now just leave the various functions as python functions.
//
// This file was auto-generated

%include "common.i"
%{
#include "twostream_interface.h"
%}

namespace FullPhysics {



class Twostream_Ls_Brdf_Supplement {

public:
  Twostream_Ls_Brdf_Supplement(const int& nbeams_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nspars_in);
  
  const blitz::Array<bool, 1>& lambertian_kernel_flag() const;
  const bool& do_shadow_effect() const;
  const bool& do_surface_emission() const;
  const int& nbeams() const;
  const int& n_user_streams() const;
  const int& n_user_relazms() const;
  const blitz::Array<double, 1>& beam_szas() const;
  const blitz::Array<double, 1>& user_angles() const;
  const blitz::Array<double, 1>& user_relazms() const;
  const double& stream_value() const;
  const int& nstreams_brdf() const;
  const int& n_brdf_kernels() const;
  const blitz::Array<int, 1>& which_brdf() const;
  const blitz::Array<double, 1>& brdf_factors() const;
  const blitz::Array<int, 1>& n_brdf_parameters() const;
  const blitz::Array<double, 2>& brdf_parameters() const;
  const int& nspars() const;
  const blitz::Array<bool, 1>& do_kernel_factor_wfs() const;
  const blitz::Array<bool, 2>& do_kernel_params_wfs() const;
  const blitz::Array<bool, 1>& do_kparams_derivs() const;
  const int& n_surface_wfs() const;
  const int& n_kernel_factor_wfs() const;
  const int& n_kernel_params_wfs() const;
  const blitz::Array<double, 2>& brdf_f_0() const;
  const blitz::Array<double, 1>& brdf_f() const;
  const blitz::Array<double, 2>& ubrdf_f() const;
  const double& emissivity() const;
  const blitz::Array<double, 3>& ls_brdf_f_0() const;
  const blitz::Array<double, 2>& ls_brdf_f() const;
  const blitz::Array<double, 3>& ls_ubrdf_f() const;
  const blitz::Array<double, 1>& ls_emissivity() const;
  
  void run();
};


class Twostream_L_Master {

public:
  Twostream_L_Master(const int& thread_in, const int& nthreads_in, const int& nlayers_in, const int& ntotal_in, const int& n_geometries_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nbeams_in, const double& earth_radius_in, const int& npars_in, const int& nspars_in);
  
  const int& thread() const;
  const bool& do_upwelling() const;
  const bool& do_dnwelling() const;
  const bool& do_plane_parallel() const;
  const bool& do_solar_sources() const;
  const bool& do_thermal_emission() const;
  const bool& do_surface_emission() const;
  const bool& do_d2s_scaling() const;
  const bool& do_brdf_surface() const;
  const bool& pure_nadir() const;
  const int& nthreads() const;
  const int& nlayers() const;
  const int& ntotal() const;
  const int& n_geometries() const;
  const double& stream_value() const;
  const int& n_user_streams() const;
  const blitz::Array<double, 1>& user_angles() const;
  const int& n_user_relazms() const;
  const blitz::Array<double, 1>& user_relazms() const;
  const double& flux_factor() const;
  const int& nbeams() const;
  const blitz::Array<double, 1>& beam_szas() const;
  const double& earth_radius() const;
  const blitz::Array<double, 1>& height_grid() const;
  const blitz::Array<double, 2>& deltau_input() const;
  const blitz::Array<double, 2>& omega_input() const;
  const blitz::Array<double, 2>& asymm_input() const;
  const blitz::Array<double, 2>& d2s_scaling() const;
  const blitz::Array<double, 1>& thermal_bb_input() const;
  const blitz::Array<double, 1>& lambertian_albedo() const;
  blitz::Array<double, 2>& brdf_f_0();
  blitz::Array<double, 1>& brdf_f();
  blitz::Array<double, 2>& ubrdf_f();
  const double& emissivity() const;
  const double& surfbb() const;
  const bool& do_sim_only() const;
  const bool& do_profile_wfs() const;
  const bool& do_column_wfs() const;
  const bool& do_surface_wfs() const;
  const int& npars() const;
  const int& nspars() const;
  const blitz::Array<bool, 1>& layer_vary_flag() const;
  const blitz::Array<int, 1>& layer_vary_number() const;
  const int& n_column_wfs() const;
  const int& n_surface_wfs() const;
  const blitz::Array<double, 3>& l_deltau_input() const;
  const blitz::Array<double, 3>& l_omega_input() const;
  const blitz::Array<double, 3>& l_asymm_input() const;
  const blitz::Array<double, 3>& l_d2s_scaling() const;
  blitz::Array<double, 3>& ls_brdf_f_0();
  blitz::Array<double, 2>& ls_brdf_f();
  blitz::Array<double, 3>& ls_ubrdf_f();
  const blitz::Array<double, 1>& ls_emiss() const;
  const blitz::Array<double, 2>& intensity_toa() const;
  const blitz::Array<double, 4>& profilewf_toa() const;
  const blitz::Array<double, 3>& columnwf_toa() const;
  const blitz::Array<double, 3>& surfacewf_toa() const;
  const blitz::Array<double, 2>& intensity_boa() const;
  const blitz::Array<double, 4>& profilewf_boa() const;
  const blitz::Array<double, 3>& columnwf_boa() const;
  const blitz::Array<double, 3>& surfacewf_boa() const;
  const int& status_inputcheck() const;
  const int& c_nmessages() const;
  std::vector< std::string > c_messages() const;
  std::vector< std::string > c_actions() const;
  const int& status_execution() const;
  std::string e_message() const;
  std::string e_trace_1() const;
  std::string e_trace_2() const;
  
  void run();
};

}
