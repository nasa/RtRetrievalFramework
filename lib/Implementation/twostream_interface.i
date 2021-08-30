// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
// This file was auto-generated

%include "common.i"

%{
#include "twostream_interface.h"
%}


  // MANUAL CHANGE
%base_import(generic_object)
  // MANUAL CHANGE

%fp_shared_ptr(FullPhysics::Twostream_Lps_Master);
%fp_shared_ptr(FullPhysics::Twostream_Ls_Brdf_Supplement);

namespace FullPhysics {



  // MANUAL CHANGE
class Twostream_Ls_Brdf_Supplement : public GenericObject {
  // MANUAL CHANGE

public:
  Twostream_Ls_Brdf_Supplement(const int& maxbeams_in, const int& max_user_streams_in, const int& max_user_obsgeoms_in, const int& maxstreams_brdf_in, const int& max_brdf_kernels_in, const int& max_brdf_parameters_in, const int& max_surfacewfs_in, const int& nbeams_in, const int& n_user_streams_in, const int& nstreams_brdf_in);
  std::string print_to_string() const;

  %python_attribute(maxbeams, int&)
  %python_attribute(max_user_streams, int&)
  %python_attribute(max_user_obsgeoms, int&)
  %python_attribute(maxstreams_brdf, int&)
  %python_attribute(max_brdf_kernels, int&)
  %python_attribute(max_brdf_parameters, int&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(do_solar_sources, bool&)
  %python_attribute(do_user_obsgeoms, bool&)
  %python_attribute(lambertian_kernel_flag, blitz::Array<bool, 1>&)
  %python_attribute(do_shadow_effect, bool&)
  %python_attribute(do_surface_emission, bool&)
  %python_attribute(nbeams, int&)
  %python_attribute(n_user_streams, int&)
  %python_attribute(n_user_obsgeoms, int&)
  %python_attribute(beam_szas, blitz::Array<double, 1>&)
  %python_attribute(user_angles, blitz::Array<double, 1>&)
  %python_attribute(user_obsgeoms, blitz::Array<double, 2>&)
  %python_attribute(stream_value, double&)
  %python_attribute(nstreams_brdf, int&)
  %python_attribute(n_brdf_kernels, int&)
  %python_attribute(which_brdf, blitz::Array<int, 1>&)
  %python_attribute(brdf_factors, blitz::Array<double, 1>&)
  %python_attribute(n_brdf_parameters, blitz::Array<int, 1>&)
  %python_attribute(brdf_parameters, blitz::Array<double, 2>&)
  %python_attribute(do_kernel_factor_wfs, blitz::Array<bool, 1>&)
  %python_attribute(do_kernel_params_wfs, blitz::Array<bool, 2>&)
  %python_attribute(do_kparams_derivs, blitz::Array<bool, 1>&)
  %python_attribute(n_surface_wfs, int&)
  %python_attribute(n_kernel_factor_wfs, int&)
  %python_attribute(n_kernel_params_wfs, int&)
  %python_attribute(brdf_f_0, blitz::Array<double, 2>&)
  %python_attribute(brdf_f, blitz::Array<double, 1>&)
  %python_attribute(ubrdf_f, blitz::Array<double, 2>&)
  %python_attribute(emissivity, double&)
  %python_attribute(ls_brdf_f_0, blitz::Array<double, 3>&)
  %python_attribute(ls_brdf_f, blitz::Array<double, 2>&)
  %python_attribute(ls_ubrdf_f, blitz::Array<double, 3>&)
  %python_attribute(ls_emissivity, blitz::Array<double, 1>&)
  %python_attribute(status_brdfsup, int&)
  
  // MANUAL CHANGE
  //%python_attribute(message, std::string&)
  //%python_attribute(action, std::string&)
  // MANUAL CHANGE
  
  void run();
  // MANUAL CHANGE
  %pickle_serialization();
  // MANUAL CHANGE
};


  // MANUAL CHANGE
class Twostream_Lps_Master : public GenericObject {
  // MANUAL CHANGE

public:
  Twostream_Lps_Master(const int& maxlayers_in, const int& maxtotal_in, const int& maxmessages_in, const int& maxbeams_in, const int& max_geometries_in, const int& max_user_streams_in, const int& max_user_relazms_in, const int& max_user_obsgeoms_in, const int& max_atmoswfs_in, const int& max_surfacewfs_in, const int& max_sleavewfs_in, const int& nlayers_in, const int& ntotal_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nbeams_in, const double& earth_radius_in, const int& n_geometries_in);
  virtual ~Twostream_Lps_Master();
  std::string print_to_string() const;

  %python_attribute(maxlayers, int&)
  %python_attribute(maxtotal, int&)
  %python_attribute(maxmessages, int&)
  %python_attribute(maxbeams, int&)
  %python_attribute(max_geometries, int&)
  %python_attribute(max_user_streams, int&)
  %python_attribute(max_user_relazms, int&)
  %python_attribute(max_user_obsgeoms, int&)
  %python_attribute(max_atmoswfs, int&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(max_sleavewfs, int&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(do_plane_parallel, bool&)
  %python_attribute(do_2s_levelout, bool&)
  %python_attribute(do_mvout_only, bool&)
  %python_attribute(do_additional_mvout, bool&)
  %python_attribute(do_solar_sources, bool&)
  %python_attribute(do_thermal_emission, bool&)
  %python_attribute(do_surface_emission, bool&)
  %python_attribute(do_d2s_scaling, bool&)
  %python_attribute(do_brdf_surface, bool&)
  %python_attribute(do_user_obsgeoms, bool&)
  %python_attribute(do_surface_leaving, bool&)
  %python_attribute(do_sl_isotropic, bool&)
  %python_attribute(do_pentadiag_inverse, bool&)
  %python_attribute(bvpindex, int&)
  %python_attribute(bvpscalefactor, double&)
  %python_attribute(taylor_order, int&)
  %python_attribute(taylor_small, double&)
  %python_attribute(tcutoff, double&)
  %python_attribute(nlayers, int&)
  %python_attribute(ntotal, int&)
  %python_attribute(stream_value, double&)
  %python_attribute(n_user_obsgeoms, int&)
  %python_attribute(user_obsgeoms, blitz::Array<double, 2>&)
  %python_attribute(n_user_streams, int&)
  %python_attribute(user_angles, blitz::Array<double, 1>&)
  %python_attribute(n_user_relazms, int&)
  %python_attribute(user_relazms, blitz::Array<double, 1>&)
  %python_attribute(flux_factor, double&)
  %python_attribute(nbeams, int&)
  %python_attribute(beam_szas, blitz::Array<double, 1>&)
  %python_attribute(earth_radius, double&)
  %python_attribute(height_grid, blitz::Array<double, 1>&)
  %python_attribute(deltau_input, blitz::Array<double, 1>&)
  %python_attribute(omega_input, blitz::Array<double, 1>&)
  %python_attribute(asymm_input, blitz::Array<double, 1>&)
  %python_attribute(d2s_scaling, blitz::Array<double, 1>&)
  %python_attribute(thermal_bb_input, blitz::Array<double, 1>&)
  %python_attribute(lambertian_albedo, double&)
  %python_attribute(brdf_f_0, blitz::Array<double, 2>&)
  %python_attribute(brdf_f, blitz::Array<double, 1>&)
  %python_attribute(ubrdf_f, blitz::Array<double, 2>&)
  %python_attribute(emissivity, double&)
  %python_attribute(surfbb, double&)
  %python_attribute(slterm_isotropic, blitz::Array<double, 1>&)
  %python_attribute(slterm_f_0, blitz::Array<double, 2>&)
  %python_attribute(do_profile_wfs, bool&)
  %python_attribute(do_surface_wfs, bool&)
  %python_attribute(do_sleave_wfs, bool&)
  %python_attribute(layer_vary_flag, blitz::Array<bool, 1>&)
  %python_attribute(layer_vary_number, blitz::Array<int, 1>&)
  %python_attribute(n_surface_wfs, int&)
  %python_attribute(n_sleave_wfs, int&)
  %python_attribute(lssl_slterm_isotropic, blitz::Array<double, 2>&)
  %python_attribute(lssl_slterm_f_0, blitz::Array<double, 3>&)
  %python_attribute(l_deltau_input, blitz::Array<double, 2>&)
  %python_attribute(l_omega_input, blitz::Array<double, 2>&)
  %python_attribute(l_asymm_input, blitz::Array<double, 2>&)
  %python_attribute(l_d2s_scaling, blitz::Array<double, 2>&)
  %python_attribute(ls_brdf_f_0, blitz::Array<double, 3>&)
  %python_attribute(ls_brdf_f, blitz::Array<double, 2>&)
  %python_attribute(ls_ubrdf_f, blitz::Array<double, 3>&)
  %python_attribute(ls_emissivity, blitz::Array<double, 1>&)
  %python_attribute(intensity_toa, blitz::Array<double, 1>&)
  %python_attribute(profilewf_toa, blitz::Array<double, 3>&)
  %python_attribute(surfacewf_toa, blitz::Array<double, 2>&)
  %python_attribute(intensity_boa, blitz::Array<double, 1>&)
  %python_attribute(profilewf_boa, blitz::Array<double, 3>&)
  %python_attribute(surfacewf_boa, blitz::Array<double, 2>&)
  %python_attribute(radlevel_up, blitz::Array<double, 2>&)
  %python_attribute(radlevel_dn, blitz::Array<double, 2>&)
  %python_attribute(n_geometries, int&)
  %python_attribute(profjaclevel_up, blitz::Array<double, 4>&)
  %python_attribute(profjaclevel_dn, blitz::Array<double, 4>&)
  %python_attribute(surfjaclevel_up, blitz::Array<double, 3>&)
  %python_attribute(surfjaclevel_dn, blitz::Array<double, 3>&)
  %python_attribute(fluxes_toa, blitz::Array<double, 2>&)
  %python_attribute(profjacfluxes_toa, blitz::Array<double, 4>&)
  %python_attribute(surfjacfluxes_toa, blitz::Array<double, 3>&)
  %python_attribute(fluxes_boa, blitz::Array<double, 2>&)
  %python_attribute(profjacfluxes_boa, blitz::Array<double, 4>&)
  %python_attribute(surfjacfluxes_boa, blitz::Array<double, 3>&)
  %python_attribute(status_inputcheck, int&)

  // MANUAL CHANGE
  const int& c_nmessages() const;
  std::vector< std::string > c_messages() const;
  std::vector< std::string > c_actions() const;
  const int& status_execution() const;
  std::string e_message() const;
  std::string e_trace_1() const;
  std::string e_trace_2() const;
  // MANUAL CHANGE

  void run();
  // MANUAL CHANGE
  %pickle_serialization();
  // MANUAL CHANGE
};

}
