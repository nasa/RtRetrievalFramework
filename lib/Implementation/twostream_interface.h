#ifndef TWOSTREAM_INTERFACE_H
#define TWOSTREAM_INTERFACE_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"


/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "twostream_ls_brdf_supplement_m" in file: "2stream_ls_brdf_supplement.F90"
//-----------------------------------------------------------------------

extern "C" {
  void twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap(const bool* lambertian_kernel_flag_in, const bool* do_shadow_effect_in, const bool* do_surface_emission_in, const int* nbeams_in, const int* n_user_streams_in, const int* n_user_relazms_in, const double* beam_szas_in, const double* user_angles_in, const double* user_relazms_in, const double* stream_value_in, const int* nstreams_brdf_in, const int* n_brdf_kernels_in, const int* which_brdf_in, const double* brdf_factors_in, const int* n_brdf_parameters_in, const double* brdf_parameters_in, const int* nspars_in, const bool* do_kernel_factor_wfs_in, const bool* do_kernel_params_wfs_in, const bool* do_kparams_derivs_in, const int* n_surface_wfs_in, const int* n_kernel_factor_wfs_in, const int* n_kernel_params_wfs_in, const double* brdf_f_0_in, const double* brdf_f_in, const double* ubrdf_f_in, const double* emissivity_in, const double* ls_brdf_f_0_in, const double* ls_brdf_f_in, const double* ls_ubrdf_f_in, const double* ls_emissivity_in);
}

class Twostream_Ls_Brdf_Supplement {

public:
  Twostream_Ls_Brdf_Supplement(const int& nbeams_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nspars_in) : nbeams_(nbeams_in), n_user_streams_(n_user_streams_in), n_user_relazms_(n_user_relazms_in), nspars_(nspars_in) 
  { 
    lambertian_kernel_flag_.reference( blitz::Array<bool, 1>(3, blitz::ColumnMajorArray<1>()) );
    lambertian_kernel_flag_ = false;
    do_shadow_effect_ = false;
    do_surface_emission_ = false;
    beam_szas_.reference( blitz::Array<double, 1>(nbeams_, blitz::ColumnMajorArray<1>()) );
    beam_szas_ = 0;
    user_angles_.reference( blitz::Array<double, 1>(n_user_streams_, blitz::ColumnMajorArray<1>()) );
    user_angles_ = 0;
    user_relazms_.reference( blitz::Array<double, 1>(n_user_relazms_, blitz::ColumnMajorArray<1>()) );
    user_relazms_ = 0;
    stream_value_ = 0;
    nstreams_brdf_ = 0;
    n_brdf_kernels_ = 0;
    which_brdf_.reference( blitz::Array<int, 1>(3, blitz::ColumnMajorArray<1>()) );
    which_brdf_ = 0;
    brdf_factors_.reference( blitz::Array<double, 1>(3, blitz::ColumnMajorArray<1>()) );
    brdf_factors_ = 0;
    n_brdf_parameters_.reference( blitz::Array<int, 1>(3, blitz::ColumnMajorArray<1>()) );
    n_brdf_parameters_ = 0;
    brdf_parameters_.reference( blitz::Array<double, 2>(3, 3, blitz::ColumnMajorArray<2>()) );
    brdf_parameters_ = 0;
    do_kernel_factor_wfs_.reference( blitz::Array<bool, 1>(3, blitz::ColumnMajorArray<1>()) );
    do_kernel_factor_wfs_ = false;
    do_kernel_params_wfs_.reference( blitz::Array<bool, 2>(3, 3, blitz::ColumnMajorArray<2>()) );
    do_kernel_params_wfs_ = false;
    do_kparams_derivs_.reference( blitz::Array<bool, 1>(3, blitz::ColumnMajorArray<1>()) );
    do_kparams_derivs_ = false;
    n_surface_wfs_ = 0;
    n_kernel_factor_wfs_ = 0;
    n_kernel_params_wfs_ = 0;
    brdf_f_0_.reference( blitz::Array<double, 2>(1-0+1, nbeams_, blitz::ColumnMajorArray<2>()) );
    brdf_f_0_ = 0;
    brdf_f_.reference( blitz::Array<double, 1>(1-0+1, blitz::ColumnMajorArray<1>()) );
    brdf_f_ = 0;
    ubrdf_f_.reference( blitz::Array<double, 2>(1-0+1, n_user_streams_, blitz::ColumnMajorArray<2>()) );
    ubrdf_f_ = 0;
    emissivity_ = 0;
    ls_brdf_f_0_.reference( blitz::Array<double, 3>(nspars_, 1-0+1, nbeams_, blitz::ColumnMajorArray<3>()) );
    ls_brdf_f_0_ = 0;
    ls_brdf_f_.reference( blitz::Array<double, 2>(nspars_, 1-0+1, blitz::ColumnMajorArray<2>()) );
    ls_brdf_f_ = 0;
    ls_ubrdf_f_.reference( blitz::Array<double, 3>(nspars_, 1-0+1, n_user_streams_, blitz::ColumnMajorArray<3>()) );
    ls_ubrdf_f_ = 0;
    ls_emissivity_.reference( blitz::Array<double, 1>(nspars_, blitz::ColumnMajorArray<1>()) );
    ls_emissivity_ = 0;
  }

  const blitz::Array<bool, 1>& lambertian_kernel_flag() const {
    return lambertian_kernel_flag_;
  }
  void lambertian_kernel_flag(const blitz::Array<bool, 1>& lambertian_kernel_flag_in) {
    lambertian_kernel_flag_ = lambertian_kernel_flag_in;
  }
  

  const bool& do_shadow_effect() const {
    return do_shadow_effect_;
  }
  void do_shadow_effect(const bool& do_shadow_effect_in) {
    do_shadow_effect_ = do_shadow_effect_in;
  }
  

  const bool& do_surface_emission() const {
    return do_surface_emission_;
  }
  void do_surface_emission(const bool& do_surface_emission_in) {
    do_surface_emission_ = do_surface_emission_in;
  }
  

  const int& nbeams() const {
    return nbeams_;
  }
  

  const int& n_user_streams() const {
    return n_user_streams_;
  }
  

  const int& n_user_relazms() const {
    return n_user_relazms_;
  }
  

  const blitz::Array<double, 1>& beam_szas() const {
    return beam_szas_;
  }
  void beam_szas(const blitz::Array<double, 1>& beam_szas_in) {
    beam_szas_ = beam_szas_in;
  }
  

  const blitz::Array<double, 1>& user_angles() const {
    return user_angles_;
  }
  void user_angles(const blitz::Array<double, 1>& user_angles_in) {
    user_angles_ = user_angles_in;
  }
  

  const blitz::Array<double, 1>& user_relazms() const {
    return user_relazms_;
  }
  void user_relazms(const blitz::Array<double, 1>& user_relazms_in) {
    user_relazms_ = user_relazms_in;
  }
  

  const double& stream_value() const {
    return stream_value_;
  }
  void stream_value(const double& stream_value_in) {
    stream_value_ = stream_value_in;
  }
  

  const int& nstreams_brdf() const {
    return nstreams_brdf_;
  }
  void nstreams_brdf(const int& nstreams_brdf_in) {
    nstreams_brdf_ = nstreams_brdf_in;
  }
  

  const int& n_brdf_kernels() const {
    return n_brdf_kernels_;
  }
  void n_brdf_kernels(const int& n_brdf_kernels_in) {
    n_brdf_kernels_ = n_brdf_kernels_in;
  }
  

  const blitz::Array<int, 1>& which_brdf() const {
    return which_brdf_;
  }
  void which_brdf(const blitz::Array<int, 1>& which_brdf_in) {
    which_brdf_ = which_brdf_in;
  }
  

  const blitz::Array<double, 1>& brdf_factors() const {
    return brdf_factors_;
  }
  void brdf_factors(const blitz::Array<double, 1>& brdf_factors_in) {
    brdf_factors_ = brdf_factors_in;
  }
  

  const blitz::Array<int, 1>& n_brdf_parameters() const {
    return n_brdf_parameters_;
  }
  void n_brdf_parameters(const blitz::Array<int, 1>& n_brdf_parameters_in) {
    n_brdf_parameters_ = n_brdf_parameters_in;
  }
  

  const blitz::Array<double, 2>& brdf_parameters() const {
    return brdf_parameters_;
  }
  void brdf_parameters(const blitz::Array<double, 2>& brdf_parameters_in) {
    brdf_parameters_ = brdf_parameters_in;
  }
  

  const int& nspars() const {
    return nspars_;
  }
  

  const blitz::Array<bool, 1>& do_kernel_factor_wfs() const {
    return do_kernel_factor_wfs_;
  }
  void do_kernel_factor_wfs(const blitz::Array<bool, 1>& do_kernel_factor_wfs_in) {
    do_kernel_factor_wfs_ = do_kernel_factor_wfs_in;
  }
  

  const blitz::Array<bool, 2>& do_kernel_params_wfs() const {
    return do_kernel_params_wfs_;
  }
  void do_kernel_params_wfs(const blitz::Array<bool, 2>& do_kernel_params_wfs_in) {
    do_kernel_params_wfs_ = do_kernel_params_wfs_in;
  }
  

  const blitz::Array<bool, 1>& do_kparams_derivs() const {
    return do_kparams_derivs_;
  }
  

  const int& n_surface_wfs() const {
    return n_surface_wfs_;
  }
  

  const int& n_kernel_factor_wfs() const {
    return n_kernel_factor_wfs_;
  }
  

  const int& n_kernel_params_wfs() const {
    return n_kernel_params_wfs_;
  }
  

  const blitz::Array<double, 2>& brdf_f_0() const {
    return brdf_f_0_;
  }
  

  const blitz::Array<double, 1>& brdf_f() const {
    return brdf_f_;
  }
  

  const blitz::Array<double, 2>& ubrdf_f() const {
    return ubrdf_f_;
  }
  

  const double& emissivity() const {
    return emissivity_;
  }
  

  const blitz::Array<double, 3>& ls_brdf_f_0() const {
    return ls_brdf_f_0_;
  }
  

  const blitz::Array<double, 2>& ls_brdf_f() const {
    return ls_brdf_f_;
  }
  

  const blitz::Array<double, 3>& ls_ubrdf_f() const {
    return ls_ubrdf_f_;
  }
  

  const blitz::Array<double, 1>& ls_emissivity() const {
    return ls_emissivity_;
  }
  

  
  void run() {
    
    
    twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap(lambertian_kernel_flag_.dataFirst(), &do_shadow_effect_, &do_surface_emission_, &nbeams_, &n_user_streams_, &n_user_relazms_, beam_szas_.dataFirst(), user_angles_.dataFirst(), user_relazms_.dataFirst(), &stream_value_, &nstreams_brdf_, &n_brdf_kernels_, which_brdf_.dataFirst(), brdf_factors_.dataFirst(), n_brdf_parameters_.dataFirst(), brdf_parameters_.dataFirst(), &nspars_, do_kernel_factor_wfs_.dataFirst(), do_kernel_params_wfs_.dataFirst(), do_kparams_derivs_.dataFirst(), &n_surface_wfs_, &n_kernel_factor_wfs_, &n_kernel_params_wfs_, brdf_f_0_.dataFirst(), brdf_f_.dataFirst(), ubrdf_f_.dataFirst(), &emissivity_, ls_brdf_f_0_.dataFirst(), ls_brdf_f_.dataFirst(), ls_ubrdf_f_.dataFirst(), ls_emissivity_.dataFirst());
    
  }

  friend std::ostream& operator<<(std::ostream &output_stream, const Twostream_Ls_Brdf_Supplement &obj) {
    output_stream << "Twostream_Ls_Brdf_Supplement:" << std::endl
      << "lambertian_kernel_flag: " << std::endl << obj.lambertian_kernel_flag()  << std::endl
      << "      do_shadow_effect: " << obj.do_shadow_effect()  << std::endl
      << "   do_surface_emission: " << obj.do_surface_emission()  << std::endl
      << "                nbeams: " << obj.nbeams()  << std::endl
      << "        n_user_streams: " << obj.n_user_streams()  << std::endl
      << "        n_user_relazms: " << obj.n_user_relazms()  << std::endl
      << "             beam_szas: " << std::endl << obj.beam_szas()  << std::endl
      << "           user_angles: " << std::endl << obj.user_angles()  << std::endl
      << "          user_relazms: " << std::endl << obj.user_relazms()  << std::endl
      << "          stream_value: " << obj.stream_value()  << std::endl
      << "         nstreams_brdf: " << obj.nstreams_brdf()  << std::endl
      << "        n_brdf_kernels: " << obj.n_brdf_kernels()  << std::endl
      << "            which_brdf: " << std::endl << obj.which_brdf()  << std::endl
      << "          brdf_factors: " << std::endl << obj.brdf_factors()  << std::endl
      << "     n_brdf_parameters: " << std::endl << obj.n_brdf_parameters()  << std::endl
      << "       brdf_parameters: " << std::endl << obj.brdf_parameters()  << std::endl
      << "                nspars: " << obj.nspars()  << std::endl
      << "  do_kernel_factor_wfs: " << std::endl << obj.do_kernel_factor_wfs()  << std::endl
      << "  do_kernel_params_wfs: " << std::endl << obj.do_kernel_params_wfs()  << std::endl
      << "     do_kparams_derivs: " << std::endl << obj.do_kparams_derivs()  << std::endl
      << "         n_surface_wfs: " << obj.n_surface_wfs()  << std::endl
      << "   n_kernel_factor_wfs: " << obj.n_kernel_factor_wfs()  << std::endl
      << "   n_kernel_params_wfs: " << obj.n_kernel_params_wfs()  << std::endl
      << "              brdf_f_0: " << std::endl << obj.brdf_f_0()  << std::endl
      << "                brdf_f: " << std::endl << obj.brdf_f()  << std::endl
      << "               ubrdf_f: " << std::endl << obj.ubrdf_f()  << std::endl
      << "            emissivity: " << obj.emissivity()  << std::endl
      << "           ls_brdf_f_0: " << std::endl << obj.ls_brdf_f_0()  << std::endl
      << "             ls_brdf_f: " << std::endl << obj.ls_brdf_f()  << std::endl
      << "            ls_ubrdf_f: " << std::endl << obj.ls_ubrdf_f()  << std::endl
      << "         ls_emissivity: " << std::endl << obj.ls_emissivity()  << std::endl;
    return output_stream;
  }

private:
  blitz::Array<bool, 1> lambertian_kernel_flag_;
  bool do_shadow_effect_;
  bool do_surface_emission_;
  int nbeams_;
  int n_user_streams_;
  int n_user_relazms_;
  blitz::Array<double, 1> beam_szas_;
  blitz::Array<double, 1> user_angles_;
  blitz::Array<double, 1> user_relazms_;
  double stream_value_;
  int nstreams_brdf_;
  int n_brdf_kernels_;
  blitz::Array<int, 1> which_brdf_;
  blitz::Array<double, 1> brdf_factors_;
  blitz::Array<int, 1> n_brdf_parameters_;
  blitz::Array<double, 2> brdf_parameters_;
  int nspars_;
  blitz::Array<bool, 1> do_kernel_factor_wfs_;
  blitz::Array<bool, 2> do_kernel_params_wfs_;
  blitz::Array<bool, 1> do_kparams_derivs_;
  int n_surface_wfs_;
  int n_kernel_factor_wfs_;
  int n_kernel_params_wfs_;
  blitz::Array<double, 2> brdf_f_0_;
  blitz::Array<double, 1> brdf_f_;
  blitz::Array<double, 2> ubrdf_f_;
  double emissivity_;
  blitz::Array<double, 3> ls_brdf_f_0_;
  blitz::Array<double, 2> ls_brdf_f_;
  blitz::Array<double, 3> ls_ubrdf_f_;
  blitz::Array<double, 1> ls_emissivity_;
};

//-----------------------------------------------------------------------
// Links to module: "twostream_l_master_m" in file: "2stream_l_master.F90"
//-----------------------------------------------------------------------

extern "C" {
  void twostream_l_master_m_twostream_l_master_wrap(const int* thread_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_plane_parallel_in, const bool* do_solar_sources_in, const bool* do_thermal_emission_in, const bool* do_surface_emission_in, const bool* do_d2s_scaling_in, const bool* do_brdf_surface_in, const bool* pure_nadir, const int* nthreads_in, const int* nlayers_in, const int* ntotal_in, const int* n_geometries_in, const double* stream_value_in, const int* n_user_streams_in, const double* user_angles_in, const int* n_user_relazms_in, const double* user_relazms_in, const double* flux_factor_in, const int* nbeams_in, const double* beam_szas_in, const double* earth_radius_in, const double* height_grid_in, const double* deltau_input_in, const double* omega_input_in, const double* asymm_input_in, const double* d2s_scaling_in, const double* thermal_bb_input_in, const double* lambertian_albedo_in, const double* brdf_f_0_in, const double* brdf_f_in, const double* ubrdf_f_in, const double* emissivity_in, const double* surfbb_in, const bool* do_sim_only_in, const bool* do_profile_wfs_in, const bool* do_column_wfs_in, const bool* do_surface_wfs_in, const int* npars_in, const int* nspars_in, const bool* layer_vary_flag_in, const int* layer_vary_number_in, const int* n_column_wfs_in, const int* n_surface_wfs_in, const double* l_deltau_input_in, const double* l_omega_input_in, const double* l_asymm_input_in, const double* l_d2s_scaling_in, const double* ls_brdf_f_0_in, const double* ls_brdf_f_in, const double* ls_ubrdf_f_in, const double* ls_emiss_in, const double* intensity_toa_in, const double* profilewf_toa_in, const double* columnwf_toa_in, const double* surfacewf_toa_in, const double* intensity_boa_in, const double* profilewf_boa_in, const double* columnwf_boa_in, const double* surfacewf_boa_in, const int* status_inputcheck_in, const int* c_nmessages_in, const int* c_messages_shape_1, const int* c_messages_len, const char* c_messages_in, const int* c_actions_shape_1, const int* c_actions_len, const char* c_actions_in, const int* status_execution_in, const int* e_message_len, const char* e_message_in, const int* e_trace_1_len, const char* e_trace_1_in, const int* e_trace_2_len, const char* e_trace_2_in);
}

class Twostream_L_Master {

public:
  Twostream_L_Master(const int& thread_in, const int& nthreads_in, const int& nlayers_in, const int& ntotal_in, const int& n_geometries_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nbeams_in, const double& earth_radius_in, const int& npars_in, const int& nspars_in) : thread_(thread_in), nthreads_(nthreads_in), nlayers_(nlayers_in), ntotal_(ntotal_in), n_geometries_(n_geometries_in), n_user_streams_(n_user_streams_in), n_user_relazms_(n_user_relazms_in), nbeams_(nbeams_in), earth_radius_(earth_radius_in), npars_(npars_in), nspars_(nspars_in) 
  { 
    do_upwelling_ = false;
    do_dnwelling_ = false;
    do_plane_parallel_ = false;
    do_solar_sources_ = false;
    do_thermal_emission_ = false;
    do_surface_emission_ = false;
    do_d2s_scaling_ = false;
    do_brdf_surface_ = false;
    pure_nadir_ = false;
    stream_value_ = 0;
    user_angles_.reference( blitz::Array<double, 1>(n_user_streams_, blitz::ColumnMajorArray<1>()) );
    user_angles_ = 0;
    user_relazms_.reference( blitz::Array<double, 1>(n_user_relazms_, blitz::ColumnMajorArray<1>()) );
    user_relazms_ = 0;
    flux_factor_ = 0;
    beam_szas_.reference( blitz::Array<double, 1>(nbeams_, blitz::ColumnMajorArray<1>()) );
    beam_szas_ = 0;
    height_grid_.reference( blitz::Array<double, 1>(nlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    height_grid_ = 0;
    deltau_input_.reference( blitz::Array<double, 2>(nlayers_, nthreads_, blitz::ColumnMajorArray<2>()) );
    deltau_input_ = 0;
    omega_input_.reference( blitz::Array<double, 2>(nlayers_, nthreads_, blitz::ColumnMajorArray<2>()) );
    omega_input_ = 0;
    asymm_input_.reference( blitz::Array<double, 2>(nlayers_, nthreads_, blitz::ColumnMajorArray<2>()) );
    asymm_input_ = 0;
    d2s_scaling_.reference( blitz::Array<double, 2>(nlayers_, nthreads_, blitz::ColumnMajorArray<2>()) );
    d2s_scaling_ = 0;
    thermal_bb_input_.reference( blitz::Array<double, 1>(nlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    thermal_bb_input_ = 0;
    lambertian_albedo_.reference( blitz::Array<double, 1>(nthreads_, blitz::ColumnMajorArray<1>()) );
    lambertian_albedo_ = 0;
    brdf_f_0_.reference( blitz::Array<double, 2>(1-0+1, nbeams_, blitz::ColumnMajorArray<2>()) );
    brdf_f_0_ = 0;
    brdf_f_.reference( blitz::Array<double, 1>(1-0+1, blitz::ColumnMajorArray<1>()) );
    brdf_f_ = 0;
    ubrdf_f_.reference( blitz::Array<double, 2>(1-0+1, n_user_streams_, blitz::ColumnMajorArray<2>()) );
    ubrdf_f_ = 0;
    emissivity_ = 0;
    surfbb_ = 0;
    do_sim_only_ = false;
    do_profile_wfs_ = false;
    do_column_wfs_ = false;
    do_surface_wfs_ = false;
    layer_vary_flag_.reference( blitz::Array<bool, 1>(nlayers_, blitz::ColumnMajorArray<1>()) );
    layer_vary_flag_ = false;
    layer_vary_number_.reference( blitz::Array<int, 1>(nlayers_, blitz::ColumnMajorArray<1>()) );
    layer_vary_number_ = 0;
    n_column_wfs_ = 0;
    n_surface_wfs_ = 0;
    l_deltau_input_.reference( blitz::Array<double, 3>(nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    l_deltau_input_ = 0;
    l_omega_input_.reference( blitz::Array<double, 3>(nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    l_omega_input_ = 0;
    l_asymm_input_.reference( blitz::Array<double, 3>(nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    l_asymm_input_ = 0;
    l_d2s_scaling_.reference( blitz::Array<double, 3>(nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    l_d2s_scaling_ = 0;
    ls_brdf_f_0_.reference( blitz::Array<double, 3>(nspars_, 1-0+1, nbeams_, blitz::ColumnMajorArray<3>()) );
    ls_brdf_f_0_ = 0;
    ls_brdf_f_.reference( blitz::Array<double, 2>(nspars_, 1-0+1, blitz::ColumnMajorArray<2>()) );
    ls_brdf_f_ = 0;
    ls_ubrdf_f_.reference( blitz::Array<double, 3>(nspars_, 1-0+1, n_user_streams_, blitz::ColumnMajorArray<3>()) );
    ls_ubrdf_f_ = 0;
    ls_emiss_.reference( blitz::Array<double, 1>(nspars_, blitz::ColumnMajorArray<1>()) );
    ls_emiss_ = 0;
    intensity_toa_.reference( blitz::Array<double, 2>(n_geometries_, nthreads_, blitz::ColumnMajorArray<2>()) );
    intensity_toa_ = 0;
    profilewf_toa_.reference( blitz::Array<double, 4>(n_geometries_, nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<4>()) );
    profilewf_toa_ = 0;
    columnwf_toa_.reference( blitz::Array<double, 3>(n_geometries_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    columnwf_toa_ = 0;
    surfacewf_toa_.reference( blitz::Array<double, 3>(n_geometries_, nspars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    surfacewf_toa_ = 0;
    intensity_boa_.reference( blitz::Array<double, 2>(n_geometries_, nthreads_, blitz::ColumnMajorArray<2>()) );
    intensity_boa_ = 0;
    profilewf_boa_.reference( blitz::Array<double, 4>(n_geometries_, nlayers_, npars_, nthreads_, blitz::ColumnMajorArray<4>()) );
    profilewf_boa_ = 0;
    columnwf_boa_.reference( blitz::Array<double, 3>(n_geometries_, npars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    columnwf_boa_ = 0;
    surfacewf_boa_.reference( blitz::Array<double, 3>(n_geometries_, nspars_, nthreads_, blitz::ColumnMajorArray<3>()) );
    surfacewf_boa_ = 0;
    status_inputcheck_ = 0;
    c_nmessages_ = 0;
    c_messages_.reference( blitz::Array<char, 2>(100, 101, blitz::ColumnMajorArray<2>()) );
    c_messages_ = '\0';
    c_actions_.reference( blitz::Array<char, 2>(100, 101, blitz::ColumnMajorArray<2>()) );
    c_actions_ = '\0';
    status_execution_ = 0;
    e_message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    e_message_ = '\0';
    e_trace_1_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    e_trace_1_ = '\0';
    e_trace_2_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    e_trace_2_ = '\0';
  }

  const int& thread() const {
    return thread_;
  }
  

  const bool& do_upwelling() const {
    return do_upwelling_;
  }
  void do_upwelling(const bool& do_upwelling_in) {
    do_upwelling_ = do_upwelling_in;
  }
  

  const bool& do_dnwelling() const {
    return do_dnwelling_;
  }
  void do_dnwelling(const bool& do_dnwelling_in) {
    do_dnwelling_ = do_dnwelling_in;
  }
  

  const bool& do_plane_parallel() const {
    return do_plane_parallel_;
  }
  void do_plane_parallel(const bool& do_plane_parallel_in) {
    do_plane_parallel_ = do_plane_parallel_in;
  }
  

  const bool& do_solar_sources() const {
    return do_solar_sources_;
  }
  void do_solar_sources(const bool& do_solar_sources_in) {
    do_solar_sources_ = do_solar_sources_in;
  }
  

  const bool& do_thermal_emission() const {
    return do_thermal_emission_;
  }
  void do_thermal_emission(const bool& do_thermal_emission_in) {
    do_thermal_emission_ = do_thermal_emission_in;
  }
  

  const bool& do_surface_emission() const {
    return do_surface_emission_;
  }
  void do_surface_emission(const bool& do_surface_emission_in) {
    do_surface_emission_ = do_surface_emission_in;
  }
  

  const bool& do_d2s_scaling() const {
    return do_d2s_scaling_;
  }
  void do_d2s_scaling(const bool& do_d2s_scaling_in) {
    do_d2s_scaling_ = do_d2s_scaling_in;
  }
  

  const bool& do_brdf_surface() const {
    return do_brdf_surface_;
  }
  void do_brdf_surface(const bool& do_brdf_surface_in) {
    do_brdf_surface_ = do_brdf_surface_in;
  }
  

  const bool& pure_nadir() const {
    return pure_nadir_;
  }
  void pure_nadir(const bool& pure_nadir_in) {
    pure_nadir_ = pure_nadir_in;
  }
  

  const int& nthreads() const {
    return nthreads_;
  }
  

  const int& nlayers() const {
    return nlayers_;
  }
  

  const int& ntotal() const {
    return ntotal_;
  }
  

  const int& n_geometries() const {
    return n_geometries_;
  }
  

  const double& stream_value() const {
    return stream_value_;
  }
  void stream_value(const double& stream_value_in) {
    stream_value_ = stream_value_in;
  }
  

  const int& n_user_streams() const {
    return n_user_streams_;
  }
  

  const blitz::Array<double, 1>& user_angles() const {
    return user_angles_;
  }
  void user_angles(const blitz::Array<double, 1>& user_angles_in) {
    user_angles_ = user_angles_in;
  }
  

  const int& n_user_relazms() const {
    return n_user_relazms_;
  }
  

  const blitz::Array<double, 1>& user_relazms() const {
    return user_relazms_;
  }
  void user_relazms(const blitz::Array<double, 1>& user_relazms_in) {
    user_relazms_ = user_relazms_in;
  }
  

  const double& flux_factor() const {
    return flux_factor_;
  }
  void flux_factor(const double& flux_factor_in) {
    flux_factor_ = flux_factor_in;
  }
  

  const int& nbeams() const {
    return nbeams_;
  }
  

  const blitz::Array<double, 1>& beam_szas() const {
    return beam_szas_;
  }
  void beam_szas(const blitz::Array<double, 1>& beam_szas_in) {
    beam_szas_ = beam_szas_in;
  }
  

  const double& earth_radius() const {
    return earth_radius_;
  }
  

  const blitz::Array<double, 1>& height_grid() const {
    return height_grid_;
  }
  void height_grid(const blitz::Array<double, 1>& height_grid_in) {
    height_grid_ = height_grid_in;
  }
  

  const blitz::Array<double, 2>& deltau_input() const {
    return deltau_input_;
  }
  void deltau_input(const blitz::Array<double, 2>& deltau_input_in) {
    deltau_input_ = deltau_input_in;
  }
  

  const blitz::Array<double, 2>& omega_input() const {
    return omega_input_;
  }
  void omega_input(const blitz::Array<double, 2>& omega_input_in) {
    omega_input_ = omega_input_in;
  }
  

  const blitz::Array<double, 2>& asymm_input() const {
    return asymm_input_;
  }
  void asymm_input(const blitz::Array<double, 2>& asymm_input_in) {
    asymm_input_ = asymm_input_in;
  }
  

  const blitz::Array<double, 2>& d2s_scaling() const {
    return d2s_scaling_;
  }
  void d2s_scaling(const blitz::Array<double, 2>& d2s_scaling_in) {
    d2s_scaling_ = d2s_scaling_in;
  }
  

  const blitz::Array<double, 1>& thermal_bb_input() const {
    return thermal_bb_input_;
  }
  void thermal_bb_input(const blitz::Array<double, 1>& thermal_bb_input_in) {
    thermal_bb_input_ = thermal_bb_input_in;
  }
  

  const blitz::Array<double, 1>& lambertian_albedo() const {
    return lambertian_albedo_;
  }
  void lambertian_albedo(const blitz::Array<double, 1>& lambertian_albedo_in) {
    lambertian_albedo_ = lambertian_albedo_in;
  }
  

  blitz::Array<double, 2>& brdf_f_0() {
    return brdf_f_0_;
  }
  const blitz::Array<double, 2>& brdf_f_0() const {
    return brdf_f_0_;
  }
  void brdf_f_0(const blitz::Array<double, 2>& brdf_f_0_in) {
    brdf_f_0_ = brdf_f_0_in;
  }
  

  blitz::Array<double, 1>& brdf_f() {
    return brdf_f_;
  }
  const blitz::Array<double, 1>& brdf_f() const {
    return brdf_f_;
  }
  void brdf_f(const blitz::Array<double, 1>& brdf_f_in) {
    brdf_f_ = brdf_f_in;
  }
  

  blitz::Array<double, 2>& ubrdf_f() {
    return ubrdf_f_;
  }
  const blitz::Array<double, 2>& ubrdf_f() const {
    return ubrdf_f_;
  }
  void ubrdf_f(const blitz::Array<double, 2>& ubrdf_f_in) {
    ubrdf_f_ = ubrdf_f_in;
  }
  

  const double& emissivity() const {
    return emissivity_;
  }
  void emissivity(const double& emissivity_in) {
    emissivity_ = emissivity_in;
  }
  

  const double& surfbb() const {
    return surfbb_;
  }
  void surfbb(const double& surfbb_in) {
    surfbb_ = surfbb_in;
  }
  

  const bool& do_sim_only() const {
    return do_sim_only_;
  }
  void do_sim_only(const bool& do_sim_only_in) {
    do_sim_only_ = do_sim_only_in;
  }
  

  const bool& do_profile_wfs() const {
    return do_profile_wfs_;
  }
  void do_profile_wfs(const bool& do_profile_wfs_in) {
    do_profile_wfs_ = do_profile_wfs_in;
  }
  

  const bool& do_column_wfs() const {
    return do_column_wfs_;
  }
  void do_column_wfs(const bool& do_column_wfs_in) {
    do_column_wfs_ = do_column_wfs_in;
  }
  

  const bool& do_surface_wfs() const {
    return do_surface_wfs_;
  }
  void do_surface_wfs(const bool& do_surface_wfs_in) {
    do_surface_wfs_ = do_surface_wfs_in;
  }
  

  const int& npars() const {
    return npars_;
  }
  

  const int& nspars() const {
    return nspars_;
  }
  

  const blitz::Array<bool, 1>& layer_vary_flag() const {
    return layer_vary_flag_;
  }
  void layer_vary_flag(const blitz::Array<bool, 1>& layer_vary_flag_in) {
    layer_vary_flag_ = layer_vary_flag_in;
  }
  

  const blitz::Array<int, 1>& layer_vary_number() const {
    return layer_vary_number_;
  }
  void layer_vary_number(const blitz::Array<int, 1>& layer_vary_number_in) {
    layer_vary_number_ = layer_vary_number_in;
  }
  

  const int& n_column_wfs() const {
    return n_column_wfs_;
  }
  void n_column_wfs(const int& n_column_wfs_in) {
    n_column_wfs_ = n_column_wfs_in;
  }
  

  const int& n_surface_wfs() const {
    return n_surface_wfs_;
  }
  void n_surface_wfs(const int& n_surface_wfs_in) {
    n_surface_wfs_ = n_surface_wfs_in;
  }
  

  const blitz::Array<double, 3>& l_deltau_input() const {
    return l_deltau_input_;
  }
  void l_deltau_input(const blitz::Array<double, 3>& l_deltau_input_in) {
    l_deltau_input_ = l_deltau_input_in;
  }
  

  const blitz::Array<double, 3>& l_omega_input() const {
    return l_omega_input_;
  }
  void l_omega_input(const blitz::Array<double, 3>& l_omega_input_in) {
    l_omega_input_ = l_omega_input_in;
  }
  

  const blitz::Array<double, 3>& l_asymm_input() const {
    return l_asymm_input_;
  }
  void l_asymm_input(const blitz::Array<double, 3>& l_asymm_input_in) {
    l_asymm_input_ = l_asymm_input_in;
  }
  

  const blitz::Array<double, 3>& l_d2s_scaling() const {
    return l_d2s_scaling_;
  }
  void l_d2s_scaling(const blitz::Array<double, 3>& l_d2s_scaling_in) {
    l_d2s_scaling_ = l_d2s_scaling_in;
  }
  

  blitz::Array<double, 3>& ls_brdf_f_0() {
    return ls_brdf_f_0_;
  }
  const blitz::Array<double, 3>& ls_brdf_f_0() const {
    return ls_brdf_f_0_;
  }
  void ls_brdf_f_0(const blitz::Array<double, 3>& ls_brdf_f_0_in) {
    ls_brdf_f_0_ = ls_brdf_f_0_in;
  }
  

  blitz::Array<double, 2>& ls_brdf_f() {
    return ls_brdf_f_;
  }
  const blitz::Array<double, 2>& ls_brdf_f() const {
    return ls_brdf_f_;
  }
  void ls_brdf_f(const blitz::Array<double, 2>& ls_brdf_f_in) {
    ls_brdf_f_ = ls_brdf_f_in;
  }
  

  blitz::Array<double, 3>& ls_ubrdf_f() {
    return ls_ubrdf_f_;
  }
  const blitz::Array<double, 3>& ls_ubrdf_f() const {
    return ls_ubrdf_f_;
  }
  void ls_ubrdf_f(const blitz::Array<double, 3>& ls_ubrdf_f_in) {
    ls_ubrdf_f_ = ls_ubrdf_f_in;
  }
  

  const blitz::Array<double, 1>& ls_emiss() const {
    return ls_emiss_;
  }
  void ls_emiss(const blitz::Array<double, 1>& ls_emiss_in) {
    ls_emiss_ = ls_emiss_in;
  }
  

  const blitz::Array<double, 2>& intensity_toa() const {
    return intensity_toa_;
  }
  void intensity_toa(const blitz::Array<double, 2>& intensity_toa_in) {
    intensity_toa_ = intensity_toa_in;
  }
  

  const blitz::Array<double, 4>& profilewf_toa() const {
    return profilewf_toa_;
  }
  void profilewf_toa(const blitz::Array<double, 4>& profilewf_toa_in) {
    profilewf_toa_ = profilewf_toa_in;
  }
  

  const blitz::Array<double, 3>& columnwf_toa() const {
    return columnwf_toa_;
  }
  void columnwf_toa(const blitz::Array<double, 3>& columnwf_toa_in) {
    columnwf_toa_ = columnwf_toa_in;
  }
  

  const blitz::Array<double, 3>& surfacewf_toa() const {
    return surfacewf_toa_;
  }
  void surfacewf_toa(const blitz::Array<double, 3>& surfacewf_toa_in) {
    surfacewf_toa_ = surfacewf_toa_in;
  }
  

  const blitz::Array<double, 2>& intensity_boa() const {
    return intensity_boa_;
  }
  void intensity_boa(const blitz::Array<double, 2>& intensity_boa_in) {
    intensity_boa_ = intensity_boa_in;
  }
  

  const blitz::Array<double, 4>& profilewf_boa() const {
    return profilewf_boa_;
  }
  void profilewf_boa(const blitz::Array<double, 4>& profilewf_boa_in) {
    profilewf_boa_ = profilewf_boa_in;
  }
  

  const blitz::Array<double, 3>& columnwf_boa() const {
    return columnwf_boa_;
  }
  void columnwf_boa(const blitz::Array<double, 3>& columnwf_boa_in) {
    columnwf_boa_ = columnwf_boa_in;
  }
  

  const blitz::Array<double, 3>& surfacewf_boa() const {
    return surfacewf_boa_;
  }
  void surfacewf_boa(const blitz::Array<double, 3>& surfacewf_boa_in) {
    surfacewf_boa_ = surfacewf_boa_in;
  }
  

  const int& status_inputcheck() const {
    return status_inputcheck_;
  }
  

  const int& c_nmessages() const {
    return c_nmessages_;
  }
  

  std::vector< std::string > c_messages() const {
    std::vector< std::string > c_messages_ret;
    for(int dim_0_idx = 0; dim_0_idx < c_messages_.extent(0); dim_0_idx++)
      c_messages_ret.push_back( std::string(std::string(c_messages_(dim_0_idx, blitz::Range::all()).begin(), c_messages_(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return c_messages_ret;
  }
  

  std::vector< std::string > c_actions() const {
    std::vector< std::string > c_actions_ret;
    for(int dim_0_idx = 0; dim_0_idx < c_actions_.extent(0); dim_0_idx++)
      c_actions_ret.push_back( std::string(std::string(c_actions_(dim_0_idx, blitz::Range::all()).begin(), c_actions_(dim_0_idx, blitz::Range::all()).end()).c_str()) );
    return c_actions_ret;
  }
  

  const int& status_execution() const {
    return status_execution_;
  }
  

  std::string e_message() const {
    std::string e_message_ret;
    e_message_ret = ( std::string(std::string(e_message_(blitz::Range::all()).begin(), e_message_(blitz::Range::all()).end()).c_str()) );
    return e_message_ret;
  }
  

  std::string e_trace_1() const {
    std::string e_trace_1_ret;
    e_trace_1_ret = ( std::string(std::string(e_trace_1_(blitz::Range::all()).begin(), e_trace_1_(blitz::Range::all()).end()).c_str()) );
    return e_trace_1_ret;
  }
  

  std::string e_trace_2() const {
    std::string e_trace_2_ret;
    e_trace_2_ret = ( std::string(std::string(e_trace_2_(blitz::Range::all()).begin(), e_trace_2_(blitz::Range::all()).end()).c_str()) );
    return e_trace_2_ret;
  }
  

  
  void run() {
    int c_messages_shape_1 = (int) c_messages_.extent(0);
    int c_messages_len = (int) c_messages_.extent(1) - 1;
    int c_actions_shape_1 = (int) c_actions_.extent(0);
    int c_actions_len = (int) c_actions_.extent(1) - 1;
    int e_message_len = (int) e_message_.extent(0) - 1;
    int e_trace_1_len = (int) e_trace_1_.extent(0) - 1;
    int e_trace_2_len = (int) e_trace_2_.extent(0) - 1;
    
    twostream_l_master_m_twostream_l_master_wrap(&thread_, &do_upwelling_, &do_dnwelling_, &do_plane_parallel_, &do_solar_sources_, &do_thermal_emission_, &do_surface_emission_, &do_d2s_scaling_, &do_brdf_surface_, &pure_nadir_, &nthreads_, &nlayers_, &ntotal_, &n_geometries_, &stream_value_, &n_user_streams_, user_angles_.dataFirst(), &n_user_relazms_, user_relazms_.dataFirst(), &flux_factor_, &nbeams_, beam_szas_.dataFirst(), &earth_radius_, height_grid_.dataFirst(), deltau_input_.dataFirst(), omega_input_.dataFirst(), asymm_input_.dataFirst(), d2s_scaling_.dataFirst(), thermal_bb_input_.dataFirst(), lambertian_albedo_.dataFirst(), brdf_f_0_.dataFirst(), brdf_f_.dataFirst(), ubrdf_f_.dataFirst(), &emissivity_, &surfbb_, &do_sim_only_, &do_profile_wfs_, &do_column_wfs_, &do_surface_wfs_, &npars_, &nspars_, layer_vary_flag_.dataFirst(), layer_vary_number_.dataFirst(), &n_column_wfs_, &n_surface_wfs_, l_deltau_input_.dataFirst(), l_omega_input_.dataFirst(), l_asymm_input_.dataFirst(), l_d2s_scaling_.dataFirst(), ls_brdf_f_0_.dataFirst(), ls_brdf_f_.dataFirst(), ls_ubrdf_f_.dataFirst(), ls_emiss_.dataFirst(), intensity_toa_.dataFirst(), profilewf_toa_.dataFirst(), columnwf_toa_.dataFirst(), surfacewf_toa_.dataFirst(), intensity_boa_.dataFirst(), profilewf_boa_.dataFirst(), columnwf_boa_.dataFirst(), surfacewf_boa_.dataFirst(), &status_inputcheck_, &c_nmessages_, &c_messages_shape_1, &c_messages_len, c_messages_.dataFirst(), &c_actions_shape_1, &c_actions_len, c_actions_.dataFirst(), &status_execution_, &e_message_len, e_message_.dataFirst(), &e_trace_1_len, e_trace_1_.dataFirst(), &e_trace_2_len, e_trace_2_.dataFirst());
    
  }
 
  friend std::ostream& operator<<(std::ostream &output_stream, const Twostream_L_Master &obj) {
    output_stream << "Twostream_L_Master:" << std::endl
      << "             thread: " << obj.thread()  << std::endl
      << "       do_upwelling: " << obj.do_upwelling()  << std::endl
      << "       do_dnwelling: " << obj.do_dnwelling()  << std::endl
      << "  do_plane_parallel: " << obj.do_plane_parallel()  << std::endl
      << "   do_solar_sources: " << obj.do_solar_sources()  << std::endl
      << "do_thermal_emission: " << obj.do_thermal_emission()  << std::endl
      << "do_surface_emission: " << obj.do_surface_emission()  << std::endl
      << "     do_d2s_scaling: " << obj.do_d2s_scaling()  << std::endl
      << "    do_brdf_surface: " << obj.do_brdf_surface()  << std::endl
      << "         pure_nadir: " << obj.pure_nadir()  << std::endl
      << "           nthreads: " << obj.nthreads()  << std::endl
      << "            nlayers: " << obj.nlayers()  << std::endl
      << "             ntotal: " << obj.ntotal()  << std::endl
      << "       n_geometries: " << obj.n_geometries()  << std::endl
      << "       stream_value: " << obj.stream_value()  << std::endl
      << "     n_user_streams: " << obj.n_user_streams()  << std::endl
      << "        user_angles: " << std::endl << obj.user_angles()  << std::endl
      << "     n_user_relazms: " << obj.n_user_relazms()  << std::endl
      << "       user_relazms: " << std::endl << obj.user_relazms()  << std::endl
      << "        flux_factor: " << obj.flux_factor()  << std::endl
      << "             nbeams: " << obj.nbeams()  << std::endl
      << "          beam_szas: " << std::endl << obj.beam_szas()  << std::endl
      << "       earth_radius: " << obj.earth_radius()  << std::endl
      << "        height_grid: " << std::endl << obj.height_grid()  << std::endl
      << "       deltau_input: " << std::endl << obj.deltau_input()  << std::endl
      << "        omega_input: " << std::endl << obj.omega_input()  << std::endl
      << "        asymm_input: " << std::endl << obj.asymm_input()  << std::endl
      << "        d2s_scaling: " << std::endl << obj.d2s_scaling()  << std::endl
      << "   thermal_bb_input: " << std::endl << obj.thermal_bb_input()  << std::endl
      << "  lambertian_albedo: " << std::endl << obj.lambertian_albedo()  << std::endl
      << "           brdf_f_0: " << std::endl << obj.brdf_f_0()  << std::endl
      << "             brdf_f: " << std::endl << obj.brdf_f()  << std::endl
      << "            ubrdf_f: " << std::endl << obj.ubrdf_f()  << std::endl
      << "         emissivity: " << obj.emissivity()  << std::endl
      << "             surfbb: " << obj.surfbb()  << std::endl
      << "        do_sim_only: " << obj.do_sim_only()  << std::endl
      << "     do_profile_wfs: " << obj.do_profile_wfs()  << std::endl
      << "      do_column_wfs: " << obj.do_column_wfs()  << std::endl
      << "     do_surface_wfs: " << obj.do_surface_wfs()  << std::endl
      << "              npars: " << obj.npars()  << std::endl
      << "             nspars: " << obj.nspars()  << std::endl
      << "    layer_vary_flag: " << std::endl << obj.layer_vary_flag()  << std::endl
      << "  layer_vary_number: " << std::endl << obj.layer_vary_number()  << std::endl
      << "       n_column_wfs: " << obj.n_column_wfs()  << std::endl
      << "      n_surface_wfs: " << obj.n_surface_wfs()  << std::endl
      << "     l_deltau_input: " << std::endl << obj.l_deltau_input()  << std::endl
      << "      l_omega_input: " << std::endl << obj.l_omega_input()  << std::endl
      << "      l_asymm_input: " << std::endl << obj.l_asymm_input()  << std::endl
      << "      l_d2s_scaling: " << std::endl << obj.l_d2s_scaling()  << std::endl
      << "        ls_brdf_f_0: " << std::endl << obj.ls_brdf_f_0()  << std::endl
      << "          ls_brdf_f: " << std::endl << obj.ls_brdf_f()  << std::endl
      << "         ls_ubrdf_f: " << std::endl << obj.ls_ubrdf_f()  << std::endl
      << "           ls_emiss: " << std::endl << obj.ls_emiss()  << std::endl
      << "      intensity_toa: " << std::endl << obj.intensity_toa()  << std::endl
      << "      profilewf_toa: " << std::endl << obj.profilewf_toa()  << std::endl
      << "       columnwf_toa: " << std::endl << obj.columnwf_toa()  << std::endl
      << "      surfacewf_toa: " << std::endl << obj.surfacewf_toa()  << std::endl
      << "      intensity_boa: " << std::endl << obj.intensity_boa()  << std::endl
      << "      profilewf_boa: " << std::endl << obj.profilewf_boa()  << std::endl
      << "       columnwf_boa: " << std::endl << obj.columnwf_boa()  << std::endl
      << "      surfacewf_boa: " << std::endl << obj.surfacewf_boa()  << std::endl
      << "  status_inputcheck: " << obj.status_inputcheck()  << std::endl
      << "        c_nmessages: " << obj.c_nmessages()  << std::endl
      << "         c_messages: " << std::endl;
    std::vector< std::string > c_messages_lcl = obj.c_messages();
    for(unsigned int idx = 0; idx < c_messages_lcl.size(); idx++)
      if ( c_messages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << c_messages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "          c_actions: " << std::endl;
    std::vector< std::string > c_actions_lcl = obj.c_actions();
    for(unsigned int idx = 0; idx < c_actions_lcl.size(); idx++)
      if ( c_actions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << c_actions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "   status_execution: " << obj.status_execution()  << std::endl
      << "          e_message: " << "\"" << obj.e_message() << "\"" << std::endl
      << "          e_trace_1: " << "\"" << obj.e_trace_1() << "\"" << std::endl
      << "          e_trace_2: " << "\"" << obj.e_trace_2() << "\"" << std::endl;
    return output_stream;
  }

private:
  int thread_;
  bool do_upwelling_;
  bool do_dnwelling_;
  bool do_plane_parallel_;
  bool do_solar_sources_;
  bool do_thermal_emission_;
  bool do_surface_emission_;
  bool do_d2s_scaling_;
  bool do_brdf_surface_;
  bool pure_nadir_;
  int nthreads_;
  int nlayers_;
  int ntotal_;
  int n_geometries_;
  double stream_value_;
  int n_user_streams_;
  blitz::Array<double, 1> user_angles_;
  int n_user_relazms_;
  blitz::Array<double, 1> user_relazms_;
  double flux_factor_;
  int nbeams_;
  blitz::Array<double, 1> beam_szas_;
  double earth_radius_;
  blitz::Array<double, 1> height_grid_;
  blitz::Array<double, 2> deltau_input_;
  blitz::Array<double, 2> omega_input_;
  blitz::Array<double, 2> asymm_input_;
  blitz::Array<double, 2> d2s_scaling_;
  blitz::Array<double, 1> thermal_bb_input_;
  blitz::Array<double, 1> lambertian_albedo_;
  blitz::Array<double, 2> brdf_f_0_;
  blitz::Array<double, 1> brdf_f_;
  blitz::Array<double, 2> ubrdf_f_;
  double emissivity_;
  double surfbb_;
  bool do_sim_only_;
  bool do_profile_wfs_;
  bool do_column_wfs_;
  bool do_surface_wfs_;
  int npars_;
  int nspars_;
  blitz::Array<bool, 1> layer_vary_flag_;
  blitz::Array<int, 1> layer_vary_number_;
  int n_column_wfs_;
  int n_surface_wfs_;
  blitz::Array<double, 3> l_deltau_input_;
  blitz::Array<double, 3> l_omega_input_;
  blitz::Array<double, 3> l_asymm_input_;
  blitz::Array<double, 3> l_d2s_scaling_;
  blitz::Array<double, 3> ls_brdf_f_0_;
  blitz::Array<double, 2> ls_brdf_f_;
  blitz::Array<double, 3> ls_ubrdf_f_;
  blitz::Array<double, 1> ls_emiss_;
  blitz::Array<double, 2> intensity_toa_;
  blitz::Array<double, 4> profilewf_toa_;
  blitz::Array<double, 3> columnwf_toa_;
  blitz::Array<double, 3> surfacewf_toa_;
  blitz::Array<double, 2> intensity_boa_;
  blitz::Array<double, 4> profilewf_boa_;
  blitz::Array<double, 3> columnwf_boa_;
  blitz::Array<double, 3> surfacewf_boa_;
  int status_inputcheck_;
  int c_nmessages_;
  blitz::Array<char, 2> c_messages_;
  blitz::Array<char, 2> c_actions_;
  int status_execution_;
  blitz::Array<char, 1> e_message_;
  blitz::Array<char, 1> e_trace_1_;
  blitz::Array<char, 1> e_trace_2_;
};



}
#endif
