#ifndef TWOSTREAM_INTERFACE_H
#define TWOSTREAM_INTERFACE_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"


/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "twostream_ls_brdf_supplement_m" in file: "2stream_ls_brdf_supplement.f90"
//-----------------------------------------------------------------------

extern "C" {
  void twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap(const int* maxbeams_in, const int* max_user_streams_in, const int* max_user_obsgeoms_in, const int* maxstreams_brdf_in, const int* max_brdf_kernels_in, const int* max_brdf_parameters_in, const int* max_surfacewfs_in, const bool* do_solar_sources_in, const bool* do_user_obsgeoms_in, const bool* lambertian_kernel_flag_in, const bool* do_shadow_effect_in, const bool* do_surface_emission_in, const int* nbeams_in, const int* n_user_streams_in, const int* n_user_obsgeoms_in, const double* beam_szas_in, const double* user_angles_in, const double* user_obsgeoms_in, const double* stream_value_in, const int* nstreams_brdf_in, const int* n_brdf_kernels_in, const int* which_brdf_in, const double* brdf_factors_in, const int* n_brdf_parameters_in, const double* brdf_parameters_in, const bool* do_kernel_factor_wfs_in, const bool* do_kernel_params_wfs_in, const bool* do_kparams_derivs_in, const int* n_surface_wfs_in, const int* n_kernel_factor_wfs_in, const int* n_kernel_params_wfs_in, const double* brdf_f_0_in, const double* brdf_f_in, const double* ubrdf_f_in, const double* emissivity_in, const double* ls_brdf_f_0_in, const double* ls_brdf_f_in, const double* ls_ubrdf_f_in, const double* ls_emissivity_in, const int* status_brdfsup_in, const int* message_len, const char* message_in, const int* action_len, const char* action_in);
}

  // MANUAL CHANGE
class Twostream_Ls_Brdf_Supplement :
    public Printable<Twostream_Ls_Brdf_Supplement> {
  // MANUAL CHANGE

public:
  Twostream_Ls_Brdf_Supplement(const int& maxbeams_in, const int& max_user_streams_in, const int& max_user_obsgeoms_in, const int& maxstreams_brdf_in, const int& max_brdf_kernels_in, const int& max_brdf_parameters_in, const int& max_surfacewfs_in, const int& nbeams_in, const int& n_user_streams_in, const int& nstreams_brdf_in) : maxbeams_(maxbeams_in), max_user_streams_(max_user_streams_in), max_user_obsgeoms_(max_user_obsgeoms_in), maxstreams_brdf_(maxstreams_brdf_in), max_brdf_kernels_(max_brdf_kernels_in), max_brdf_parameters_(max_brdf_parameters_in), max_surfacewfs_(max_surfacewfs_in), nbeams_(nbeams_in), n_user_streams_(n_user_streams_in), nstreams_brdf_(nstreams_brdf_in) 
  { 
    do_solar_sources_ = false;
    do_user_obsgeoms_ = false;
    lambertian_kernel_flag_.reference( blitz::Array<bool, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    lambertian_kernel_flag_ = false;
    do_shadow_effect_ = false;
    do_surface_emission_ = false;
    n_user_obsgeoms_ = 0;
    beam_szas_.reference( blitz::Array<double, 1>(maxbeams_, blitz::ColumnMajorArray<1>()) );
    beam_szas_ = 0;
    user_angles_.reference( blitz::Array<double, 1>(max_user_streams_, blitz::ColumnMajorArray<1>()) );
    user_angles_ = 0;
    user_obsgeoms_.reference( blitz::Array<double, 2>(max_user_obsgeoms_, 3, blitz::ColumnMajorArray<2>()) );
    user_obsgeoms_ = 0;
    stream_value_ = 0;
    n_brdf_kernels_ = 0;
    which_brdf_.reference( blitz::Array<int, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    which_brdf_ = 0;
    brdf_factors_.reference( blitz::Array<double, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    brdf_factors_ = 0;
    n_brdf_parameters_.reference( blitz::Array<int, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    n_brdf_parameters_ = 0;
    brdf_parameters_.reference( blitz::Array<double, 2>(max_brdf_kernels_, max_brdf_parameters_, blitz::ColumnMajorArray<2>()) );
    brdf_parameters_ = 0;
    do_kernel_factor_wfs_.reference( blitz::Array<bool, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    do_kernel_factor_wfs_ = false;
    do_kernel_params_wfs_.reference( blitz::Array<bool, 2>(max_brdf_kernels_, max_brdf_parameters_, blitz::ColumnMajorArray<2>()) );
    do_kernel_params_wfs_ = false;
    do_kparams_derivs_.reference( blitz::Array<bool, 1>(max_brdf_kernels_, blitz::ColumnMajorArray<1>()) );
    do_kparams_derivs_ = false;
    n_surface_wfs_ = 0;
    n_kernel_factor_wfs_ = 0;
    n_kernel_params_wfs_ = 0;
    brdf_f_0_.reference( blitz::Array<double, 2>(1-0+1, maxbeams_, blitz::ColumnMajorArray<2>()) );
    brdf_f_0_ = 0;
    brdf_f_.reference( blitz::Array<double, 1>(1-0+1, blitz::ColumnMajorArray<1>()) );
    brdf_f_ = 0;
    ubrdf_f_.reference( blitz::Array<double, 2>(1-0+1, max_user_streams_, blitz::ColumnMajorArray<2>()) );
    ubrdf_f_ = 0;
    emissivity_ = 0;
    ls_brdf_f_0_.reference( blitz::Array<double, 3>(max_surfacewfs_, 1-0+1, maxbeams_, blitz::ColumnMajorArray<3>()) );
    ls_brdf_f_0_ = 0;
    ls_brdf_f_.reference( blitz::Array<double, 2>(max_surfacewfs_, 1-0+1, blitz::ColumnMajorArray<2>()) );
    ls_brdf_f_ = 0;
    ls_ubrdf_f_.reference( blitz::Array<double, 3>(max_surfacewfs_, 1-0+1, max_user_streams_, blitz::ColumnMajorArray<3>()) );
    ls_ubrdf_f_ = 0;
    ls_emissivity_.reference( blitz::Array<double, 1>(max_surfacewfs_, blitz::ColumnMajorArray<1>()) );
    ls_emissivity_ = 0;
    status_brdfsup_ = 0;
    message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    message_ = '\0';
    action_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    action_ = '\0';
    // Initialize type pointers
    
  }

  virtual ~Twostream_Ls_Brdf_Supplement() = default;

  const int& maxbeams() const {
    return maxbeams_;
  }

  

  const int& max_user_streams() const {
    return max_user_streams_;
  }

  

  const int& max_user_obsgeoms() const {
    return max_user_obsgeoms_;
  }

  

  const int& maxstreams_brdf() const {
    return maxstreams_brdf_;
  }

  

  const int& max_brdf_kernels() const {
    return max_brdf_kernels_;
  }

  

  const int& max_brdf_parameters() const {
    return max_brdf_parameters_;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const bool& do_solar_sources() const {
    return do_solar_sources_;
  }

  void do_solar_sources(const bool& do_solar_sources_in) {
    do_solar_sources_ = do_solar_sources_in;
  }

  

  const bool& do_user_obsgeoms() const {
    return do_user_obsgeoms_;
  }

  void do_user_obsgeoms(const bool& do_user_obsgeoms_in) {
    do_user_obsgeoms_ = do_user_obsgeoms_in;
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

  

  const int& n_user_obsgeoms() const {
    return n_user_obsgeoms_;
  }

  void n_user_obsgeoms(const int& n_user_obsgeoms_in) {
    n_user_obsgeoms_ = n_user_obsgeoms_in;
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

  

  const blitz::Array<double, 2>& user_obsgeoms() const {
    return user_obsgeoms_;
  }

  void user_obsgeoms(const blitz::Array<double, 2>& user_obsgeoms_in) {
    user_obsgeoms_ = user_obsgeoms_in;
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

  

  const int& status_brdfsup() const {
    return status_brdfsup_;
  }

  

  std::string message() const {
    std::string message_ret;
    message_ret = ( std::string(std::string(message_(blitz::Range::all()).begin(), message_(blitz::Range::all()).end()).c_str()) );
    return message_ret;
  }

  

  std::string action() const {
    std::string action_ret;
    action_ret = ( std::string(std::string(action_(blitz::Range::all()).begin(), action_(blitz::Range::all()).end()).c_str()) );
    return action_ret;
  }

  

  
  void run() {
    int message_len = (int) message_.extent(0) - 1;
    int action_len = (int) action_.extent(0) - 1;
    
    twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap(&maxbeams_, &max_user_streams_, &max_user_obsgeoms_, &maxstreams_brdf_, &max_brdf_kernels_, &max_brdf_parameters_, &max_surfacewfs_, &do_solar_sources_, &do_user_obsgeoms_, lambertian_kernel_flag_.dataFirst(), &do_shadow_effect_, &do_surface_emission_, &nbeams_, &n_user_streams_, &n_user_obsgeoms_, beam_szas_.dataFirst(), user_angles_.dataFirst(), user_obsgeoms_.dataFirst(), &stream_value_, &nstreams_brdf_, &n_brdf_kernels_, which_brdf_.dataFirst(), brdf_factors_.dataFirst(), n_brdf_parameters_.dataFirst(), brdf_parameters_.dataFirst(), do_kernel_factor_wfs_.dataFirst(), do_kernel_params_wfs_.dataFirst(), do_kparams_derivs_.dataFirst(), &n_surface_wfs_, &n_kernel_factor_wfs_, &n_kernel_params_wfs_, brdf_f_0_.dataFirst(), brdf_f_.dataFirst(), ubrdf_f_.dataFirst(), &emissivity_, ls_brdf_f_0_.dataFirst(), ls_brdf_f_.dataFirst(), ls_ubrdf_f_.dataFirst(), ls_emissivity_.dataFirst(), &status_brdfsup_, &message_len, message_.dataFirst(), &action_len, action_.dataFirst());
    
  }

  // MANUAL CHANGE
  void print(std::ostream &output_stream) const {
    output_stream << "Twostream_Ls_Brdf_Supplement:" << std::endl
      << "              maxbeams: " << maxbeams()  << std::endl
      << "      max_user_streams: " << max_user_streams()  << std::endl
      << "     max_user_obsgeoms: " << max_user_obsgeoms()  << std::endl
      << "       maxstreams_brdf: " << maxstreams_brdf()  << std::endl
      << "      max_brdf_kernels: " << max_brdf_kernels()  << std::endl
      << "   max_brdf_parameters: " << max_brdf_parameters()  << std::endl
      << "        max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "      do_solar_sources: " << do_solar_sources()  << std::endl
      << "      do_user_obsgeoms: " << do_user_obsgeoms()  << std::endl
      << "lambertian_kernel_flag: " << std::endl << lambertian_kernel_flag()  << std::endl
      << "      do_shadow_effect: " << do_shadow_effect()  << std::endl
      << "   do_surface_emission: " << do_surface_emission()  << std::endl
      << "                nbeams: " << nbeams()  << std::endl
      << "        n_user_streams: " << n_user_streams()  << std::endl
      << "       n_user_obsgeoms: " << n_user_obsgeoms()  << std::endl
      << "             beam_szas: " << std::endl << beam_szas()  << std::endl
      << "           user_angles: " << std::endl << user_angles()  << std::endl
      << "         user_obsgeoms: " << std::endl << user_obsgeoms()  << std::endl
      << "          stream_value: " << stream_value()  << std::endl
      << "         nstreams_brdf: " << nstreams_brdf()  << std::endl
      << "        n_brdf_kernels: " << n_brdf_kernels()  << std::endl
      << "            which_brdf: " << std::endl << which_brdf()  << std::endl
      << "          brdf_factors: " << std::endl << brdf_factors()  << std::endl
      << "     n_brdf_parameters: " << std::endl << n_brdf_parameters()  << std::endl
      << "       brdf_parameters: " << std::endl << brdf_parameters()  << std::endl
      << "  do_kernel_factor_wfs: " << std::endl << do_kernel_factor_wfs()  << std::endl
      << "  do_kernel_params_wfs: " << std::endl << do_kernel_params_wfs()  << std::endl
      << "     do_kparams_derivs: " << std::endl << do_kparams_derivs()  << std::endl
      << "         n_surface_wfs: " << n_surface_wfs()  << std::endl
      << "   n_kernel_factor_wfs: " << n_kernel_factor_wfs()  << std::endl
      << "   n_kernel_params_wfs: " << n_kernel_params_wfs()  << std::endl
      << "              brdf_f_0: " << std::endl << brdf_f_0()  << std::endl
      << "                brdf_f: " << std::endl << brdf_f()  << std::endl
      << "               ubrdf_f: " << std::endl << ubrdf_f()  << std::endl
      << "            emissivity: " << emissivity()  << std::endl
      << "           ls_brdf_f_0: " << std::endl << ls_brdf_f_0()  << std::endl
      << "             ls_brdf_f: " << std::endl << ls_brdf_f()  << std::endl
      << "            ls_ubrdf_f: " << std::endl << ls_ubrdf_f()  << std::endl
      << "         ls_emissivity: " << std::endl << ls_emissivity()  << std::endl
      << "        status_brdfsup: " << status_brdfsup()  << std::endl
      << "               message: " << "\"" << message() << "\"" << std::endl
      << "                action: " << "\"" << action() << "\"" << std::endl;

  }
  // MANUAL CHANGE

private:
  int maxbeams_;
  int max_user_streams_;
  int max_user_obsgeoms_;
  int maxstreams_brdf_;
  int max_brdf_kernels_;
  int max_brdf_parameters_;
  int max_surfacewfs_;
  bool do_solar_sources_;
  bool do_user_obsgeoms_;
  blitz::Array<bool, 1> lambertian_kernel_flag_;
  bool do_shadow_effect_;
  bool do_surface_emission_;
  int nbeams_;
  int n_user_streams_;
  int n_user_obsgeoms_;
  blitz::Array<double, 1> beam_szas_;
  blitz::Array<double, 1> user_angles_;
  blitz::Array<double, 2> user_obsgeoms_;
  double stream_value_;
  int nstreams_brdf_;
  int n_brdf_kernels_;
  blitz::Array<int, 1> which_brdf_;
  blitz::Array<double, 1> brdf_factors_;
  blitz::Array<int, 1> n_brdf_parameters_;
  blitz::Array<double, 2> brdf_parameters_;
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
  int status_brdfsup_;
  blitz::Array<char, 1> message_;
  blitz::Array<char, 1> action_;
  // MANUAL CHANGE
  Twostream_Ls_Brdf_Supplement() {}
  // MANUAL CHANGE
};

//-----------------------------------------------------------------------
// Links to module: "twostream_lps_master_m" in file: "2stream_lps_master.f90"
//-----------------------------------------------------------------------

extern "C" {
  void twostream_lps_master_m_twostream_lps_master_wrap(const int* maxlayers_in, const int* maxtotal_in, const int* maxmessages_in, const int* maxbeams_in, const int* max_geometries_in, const int* max_user_streams_in, const int* max_user_relazms_in, const int* max_user_obsgeoms_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const int* max_sleavewfs_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_plane_parallel_in, const bool* do_2s_levelout_in, const bool* do_mvout_only_in, const bool* do_additional_mvout_in, const bool* do_solar_sources_in, const bool* do_thermal_emission_in, const bool* do_surface_emission_in, const bool* do_d2s_scaling_in, const bool* do_brdf_surface_in, const bool* do_user_obsgeoms_in, const bool* do_surface_leaving_in, const bool* do_sl_isotropic_in, const bool* do_pentadiag_inverse_in, const int* bvpindex_in, const double* bvpscalefactor_in, const int* taylor_order_in, const double* taylor_small_in, const double* tcutoff_in, const int* nlayers_in, const int* ntotal_in, const double* stream_value_in, const int* n_user_obsgeoms_in, const double* user_obsgeoms_in, const int* n_user_streams_in, const double* user_angles_in, const int* n_user_relazms_in, const double* user_relazms_in, const double* flux_factor_in, const int* nbeams_in, const double* beam_szas_in, const double* earth_radius_in, const double* height_grid_in, const double* deltau_input_in, const double* omega_input_in, const double* asymm_input_in, const double* d2s_scaling_in, const double* thermal_bb_input_in, const double* lambertian_albedo_in, const double* brdf_f_0_in, const double* brdf_f_in, const double* ubrdf_f_in, const double* emissivity_in, const double* surfbb_in, const double* slterm_isotropic_in, const double* slterm_f_0_in, const bool* do_profile_wfs_in, const bool* do_surface_wfs_in, const bool* do_sleave_wfs_in, const bool* layer_vary_flag_in, const int* layer_vary_number_in, const int* n_surface_wfs_in, const int* n_sleave_wfs_in, const double* lssl_slterm_isotropic_in, const double* lssl_slterm_f_0_in, const double* l_deltau_input_in, const double* l_omega_input_in, const double* l_asymm_input_in, const double* l_d2s_scaling_in, const double* ls_brdf_f_0_in, const double* ls_brdf_f_in, const double* ls_ubrdf_f_in, const double* ls_emissivity_in, const double* intensity_toa_in, const double* profilewf_toa_in, const double* surfacewf_toa_in, const double* intensity_boa_in, const double* profilewf_boa_in, const double* surfacewf_boa_in, const double* radlevel_up_in, const double* radlevel_dn_in, const int* n_geometries_in, const double* profjaclevel_up_in, const double* profjaclevel_dn_in, const double* surfjaclevel_up_in, const double* surfjaclevel_dn_in, const double* fluxes_toa_in, const double* profjacfluxes_toa_in, const double* surfjacfluxes_toa_in, const double* fluxes_boa_in, const double* profjacfluxes_boa_in, const double* surfjacfluxes_boa_in, const int* status_inputcheck_in, const int* c_nmessages_in, const int* c_messages_shape_1, const int* c_messages_len, const char* c_messages_in, const int* c_actions_shape_1, const int* c_actions_len, const char* c_actions_in, const int* status_execution_in, const int* e_message_len, const char* e_message_in, const int* e_trace_1_len, const char* e_trace_1_in, const int* e_trace_2_len, const char* e_trace_2_in);
}

  // MANUAL CHANGE
class Twostream_Lps_Master : public Printable<Twostream_Lps_Master> {
  // MANUAL CHANGE

public:
  Twostream_Lps_Master(const int& maxlayers_in, const int& maxtotal_in, const int& maxmessages_in, const int& maxbeams_in, const int& max_geometries_in, const int& max_user_streams_in, const int& max_user_relazms_in, const int& max_user_obsgeoms_in, const int& max_atmoswfs_in, const int& max_surfacewfs_in, const int& max_sleavewfs_in, const int& nlayers_in, const int& ntotal_in, const int& n_user_streams_in, const int& n_user_relazms_in, const int& nbeams_in, const double& earth_radius_in, const int& n_geometries_in) : maxlayers_(maxlayers_in), maxtotal_(maxtotal_in), maxmessages_(maxmessages_in), maxbeams_(maxbeams_in), max_geometries_(max_geometries_in), max_user_streams_(max_user_streams_in), max_user_relazms_(max_user_relazms_in), max_user_obsgeoms_(max_user_obsgeoms_in), max_atmoswfs_(max_atmoswfs_in), max_surfacewfs_(max_surfacewfs_in), max_sleavewfs_(max_sleavewfs_in), nlayers_(nlayers_in), ntotal_(ntotal_in), n_user_streams_(n_user_streams_in), n_user_relazms_(n_user_relazms_in), nbeams_(nbeams_in), earth_radius_(earth_radius_in), n_geometries_(n_geometries_in) 
  { 
    do_upwelling_ = false;
    do_dnwelling_ = false;
    do_plane_parallel_ = false;
    do_2s_levelout_ = false;
    do_mvout_only_ = false;
    do_additional_mvout_ = false;
    do_solar_sources_ = false;
    do_thermal_emission_ = false;
    do_surface_emission_ = false;
    do_d2s_scaling_ = false;
    do_brdf_surface_ = false;
    do_user_obsgeoms_ = false;
    do_surface_leaving_ = false;
    do_sl_isotropic_ = false;
    do_pentadiag_inverse_ = false;
    bvpindex_ = 0;
    bvpscalefactor_ = 0;
    taylor_order_ = 0;
    taylor_small_ = 0;
    tcutoff_ = 0;
    stream_value_ = 0;
    n_user_obsgeoms_ = 0;
    user_obsgeoms_.reference( blitz::Array<double, 2>(max_user_obsgeoms_, 3, blitz::ColumnMajorArray<2>()) );
    user_obsgeoms_ = 0;
    user_angles_.reference( blitz::Array<double, 1>(max_user_streams_, blitz::ColumnMajorArray<1>()) );
    user_angles_ = 0;
    user_relazms_.reference( blitz::Array<double, 1>(max_user_relazms_, blitz::ColumnMajorArray<1>()) );
    user_relazms_ = 0;
    flux_factor_ = 0;
    beam_szas_.reference( blitz::Array<double, 1>(maxbeams_, blitz::ColumnMajorArray<1>()) );
    beam_szas_ = 0;
    height_grid_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    height_grid_ = 0;
    deltau_input_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltau_input_ = 0;
    omega_input_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    omega_input_ = 0;
    asymm_input_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    asymm_input_ = 0;
    d2s_scaling_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    d2s_scaling_ = 0;
    thermal_bb_input_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    thermal_bb_input_ = 0;
    lambertian_albedo_ = 0;
    brdf_f_0_.reference( blitz::Array<double, 2>(1-0+1, maxbeams_, blitz::ColumnMajorArray<2>()) );
    brdf_f_0_ = 0;
    brdf_f_.reference( blitz::Array<double, 1>(1-0+1, blitz::ColumnMajorArray<1>()) );
    brdf_f_ = 0;
    ubrdf_f_.reference( blitz::Array<double, 2>(1-0+1, max_user_streams_, blitz::ColumnMajorArray<2>()) );
    ubrdf_f_ = 0;
    emissivity_ = 0;
    surfbb_ = 0;
    slterm_isotropic_.reference( blitz::Array<double, 1>(maxbeams_, blitz::ColumnMajorArray<1>()) );
    slterm_isotropic_ = 0;
    slterm_f_0_.reference( blitz::Array<double, 2>(1-0+1, maxbeams_, blitz::ColumnMajorArray<2>()) );
    slterm_f_0_ = 0;
    do_profile_wfs_ = false;
    do_surface_wfs_ = false;
    do_sleave_wfs_ = false;
    layer_vary_flag_.reference( blitz::Array<bool, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    layer_vary_flag_ = false;
    layer_vary_number_.reference( blitz::Array<int, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    layer_vary_number_ = 0;
    n_surface_wfs_ = 0;
    n_sleave_wfs_ = 0;
    lssl_slterm_isotropic_.reference( blitz::Array<double, 2>(max_sleavewfs_, maxbeams_, blitz::ColumnMajorArray<2>()) );
    lssl_slterm_isotropic_ = 0;
    lssl_slterm_f_0_.reference( blitz::Array<double, 3>(max_sleavewfs_, 1-0+1, maxbeams_, blitz::ColumnMajorArray<3>()) );
    lssl_slterm_f_0_ = 0;
    l_deltau_input_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_deltau_input_ = 0;
    l_omega_input_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_omega_input_ = 0;
    l_asymm_input_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_asymm_input_ = 0;
    l_d2s_scaling_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_d2s_scaling_ = 0;
    ls_brdf_f_0_.reference( blitz::Array<double, 3>(max_surfacewfs_, 1-0+1, maxbeams_, blitz::ColumnMajorArray<3>()) );
    ls_brdf_f_0_ = 0;
    ls_brdf_f_.reference( blitz::Array<double, 2>(max_surfacewfs_, 1-0+1, blitz::ColumnMajorArray<2>()) );
    ls_brdf_f_ = 0;
    ls_ubrdf_f_.reference( blitz::Array<double, 3>(max_surfacewfs_, 1-0+1, max_user_streams_, blitz::ColumnMajorArray<3>()) );
    ls_ubrdf_f_ = 0;
    ls_emissivity_.reference( blitz::Array<double, 1>(max_surfacewfs_, blitz::ColumnMajorArray<1>()) );
    ls_emissivity_ = 0;
    intensity_toa_.reference( blitz::Array<double, 1>(max_geometries_, blitz::ColumnMajorArray<1>()) );
    intensity_toa_ = 0;
    profilewf_toa_.reference( blitz::Array<double, 3>(max_geometries_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    profilewf_toa_ = 0;
    surfacewf_toa_.reference( blitz::Array<double, 2>(max_geometries_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    surfacewf_toa_ = 0;
    intensity_boa_.reference( blitz::Array<double, 1>(max_geometries_, blitz::ColumnMajorArray<1>()) );
    intensity_boa_ = 0;
    profilewf_boa_.reference( blitz::Array<double, 3>(max_geometries_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    profilewf_boa_ = 0;
    surfacewf_boa_.reference( blitz::Array<double, 2>(max_geometries_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    surfacewf_boa_ = 0;
    radlevel_up_.reference( blitz::Array<double, 2>(max_geometries_, maxlayers_-0+1, blitz::ColumnMajorArray<2>()) );
    radlevel_up_ = 0;
    radlevel_dn_.reference( blitz::Array<double, 2>(max_geometries_, maxlayers_-0+1, blitz::ColumnMajorArray<2>()) );
    radlevel_dn_ = 0;
    profjaclevel_up_.reference( blitz::Array<double, 4>(max_geometries_, maxlayers_-0+1, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    profjaclevel_up_ = 0;
    profjaclevel_dn_.reference( blitz::Array<double, 4>(max_geometries_, maxlayers_-0+1, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    profjaclevel_dn_ = 0;
    surfjaclevel_up_.reference( blitz::Array<double, 3>(max_geometries_, maxlayers_-0+1, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    surfjaclevel_up_ = 0;
    surfjaclevel_dn_.reference( blitz::Array<double, 3>(max_geometries_, maxlayers_-0+1, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    surfjaclevel_dn_ = 0;
    fluxes_toa_.reference( blitz::Array<double, 2>(maxbeams_, 2, blitz::ColumnMajorArray<2>()) );
    fluxes_toa_ = 0;
    profjacfluxes_toa_.reference( blitz::Array<double, 4>(maxbeams_, 2, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    profjacfluxes_toa_ = 0;
    surfjacfluxes_toa_.reference( blitz::Array<double, 3>(maxbeams_, 2, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    surfjacfluxes_toa_ = 0;
    fluxes_boa_.reference( blitz::Array<double, 2>(maxbeams_, 2, blitz::ColumnMajorArray<2>()) );
    fluxes_boa_ = 0;
    profjacfluxes_boa_.reference( blitz::Array<double, 4>(maxbeams_, 2, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    profjacfluxes_boa_ = 0;
    surfjacfluxes_boa_.reference( blitz::Array<double, 3>(maxbeams_, 2, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    surfjacfluxes_boa_ = 0;
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
    // Initialize type pointers
    
  }

  virtual ~Twostream_Lps_Master() = default;

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxtotal() const {
    return maxtotal_;
  }

  

  const int& maxmessages() const {
    return maxmessages_;
  }

  

  const int& maxbeams() const {
    return maxbeams_;
  }

  

  const int& max_geometries() const {
    return max_geometries_;
  }

  

  const int& max_user_streams() const {
    return max_user_streams_;
  }

  

  const int& max_user_relazms() const {
    return max_user_relazms_;
  }

  

  const int& max_user_obsgeoms() const {
    return max_user_obsgeoms_;
  }

  

  const int& max_atmoswfs() const {
    return max_atmoswfs_;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const int& max_sleavewfs() const {
    return max_sleavewfs_;
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

  

  const bool& do_2s_levelout() const {
    return do_2s_levelout_;
  }

  void do_2s_levelout(const bool& do_2s_levelout_in) {
    do_2s_levelout_ = do_2s_levelout_in;
  }

  

  const bool& do_mvout_only() const {
    return do_mvout_only_;
  }

  void do_mvout_only(const bool& do_mvout_only_in) {
    do_mvout_only_ = do_mvout_only_in;
  }

  

  const bool& do_additional_mvout() const {
    return do_additional_mvout_;
  }

  void do_additional_mvout(const bool& do_additional_mvout_in) {
    do_additional_mvout_ = do_additional_mvout_in;
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

  

  const bool& do_user_obsgeoms() const {
    return do_user_obsgeoms_;
  }

  void do_user_obsgeoms(const bool& do_user_obsgeoms_in) {
    do_user_obsgeoms_ = do_user_obsgeoms_in;
  }

  

  const bool& do_surface_leaving() const {
    return do_surface_leaving_;
  }

  void do_surface_leaving(const bool& do_surface_leaving_in) {
    do_surface_leaving_ = do_surface_leaving_in;
  }

  

  const bool& do_sl_isotropic() const {
    return do_sl_isotropic_;
  }

  void do_sl_isotropic(const bool& do_sl_isotropic_in) {
    do_sl_isotropic_ = do_sl_isotropic_in;
  }

  

  const bool& do_pentadiag_inverse() const {
    return do_pentadiag_inverse_;
  }

  void do_pentadiag_inverse(const bool& do_pentadiag_inverse_in) {
    do_pentadiag_inverse_ = do_pentadiag_inverse_in;
  }

  

  const int& bvpindex() const {
    return bvpindex_;
  }

  void bvpindex(const int& bvpindex_in) {
    bvpindex_ = bvpindex_in;
  }

  

  const double& bvpscalefactor() const {
    return bvpscalefactor_;
  }

  void bvpscalefactor(const double& bvpscalefactor_in) {
    bvpscalefactor_ = bvpscalefactor_in;
  }

  

  const int& taylor_order() const {
    return taylor_order_;
  }

  void taylor_order(const int& taylor_order_in) {
    taylor_order_ = taylor_order_in;
  }

  

  const double& taylor_small() const {
    return taylor_small_;
  }

  void taylor_small(const double& taylor_small_in) {
    taylor_small_ = taylor_small_in;
  }

  

  const double& tcutoff() const {
    return tcutoff_;
  }

  void tcutoff(const double& tcutoff_in) {
    tcutoff_ = tcutoff_in;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const int& ntotal() const {
    return ntotal_;
  }

  

  const double& stream_value() const {
    return stream_value_;
  }

  void stream_value(const double& stream_value_in) {
    stream_value_ = stream_value_in;
  }

  

  const int& n_user_obsgeoms() const {
    return n_user_obsgeoms_;
  }

  void n_user_obsgeoms(const int& n_user_obsgeoms_in) {
    n_user_obsgeoms_ = n_user_obsgeoms_in;
  }

  

  const blitz::Array<double, 2>& user_obsgeoms() const {
    return user_obsgeoms_;
  }

  void user_obsgeoms(const blitz::Array<double, 2>& user_obsgeoms_in) {
    user_obsgeoms_ = user_obsgeoms_in;
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

  

  const blitz::Array<double, 1>& deltau_input() const {
    return deltau_input_;
  }

  void deltau_input(const blitz::Array<double, 1>& deltau_input_in) {
    deltau_input_ = deltau_input_in;
  }

  

  const blitz::Array<double, 1>& omega_input() const {
    return omega_input_;
  }

  void omega_input(const blitz::Array<double, 1>& omega_input_in) {
    omega_input_ = omega_input_in;
  }

  

  const blitz::Array<double, 1>& asymm_input() const {
    return asymm_input_;
  }

  void asymm_input(const blitz::Array<double, 1>& asymm_input_in) {
    asymm_input_ = asymm_input_in;
  }

  

  const blitz::Array<double, 1>& d2s_scaling() const {
    return d2s_scaling_;
  }

  void d2s_scaling(const blitz::Array<double, 1>& d2s_scaling_in) {
    d2s_scaling_ = d2s_scaling_in;
  }

  

  const blitz::Array<double, 1>& thermal_bb_input() const {
    return thermal_bb_input_;
  }

  void thermal_bb_input(const blitz::Array<double, 1>& thermal_bb_input_in) {
    thermal_bb_input_ = thermal_bb_input_in;
  }

  

  const double& lambertian_albedo() const {
    return lambertian_albedo_;
  }

  void lambertian_albedo(const double& lambertian_albedo_in) {
    lambertian_albedo_ = lambertian_albedo_in;
  }

  
  // MANUAL EDIT //
  blitz::Array<double, 2>& brdf_f_0() {
    return brdf_f_0_;
  }
  // MANUAL EDIT //


  const blitz::Array<double, 2>& brdf_f_0() const {
    return brdf_f_0_;
  }

  void brdf_f_0(const blitz::Array<double, 2>& brdf_f_0_in) {
    brdf_f_0_ = brdf_f_0_in;
  }

  
  // MANUAL EDIT //
  blitz::Array<double, 1>& brdf_f() {
    return brdf_f_;
  }
  // MANUAL EDIT //

  const blitz::Array<double, 1>& brdf_f() const {
    return brdf_f_;
  }

  void brdf_f(const blitz::Array<double, 1>& brdf_f_in) {
    brdf_f_ = brdf_f_in;
  }

  // MANUAL EDIT //
  blitz::Array<double, 2>& ubrdf_f() {
    return ubrdf_f_;
  }
  // MANUAL EDIT //

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

  

  const blitz::Array<double, 1>& slterm_isotropic() const {
    return slterm_isotropic_;
  }

  void slterm_isotropic(const blitz::Array<double, 1>& slterm_isotropic_in) {
    slterm_isotropic_ = slterm_isotropic_in;
  }

  

  const blitz::Array<double, 2>& slterm_f_0() const {
    return slterm_f_0_;
  }

  void slterm_f_0(const blitz::Array<double, 2>& slterm_f_0_in) {
    slterm_f_0_ = slterm_f_0_in;
  }

  

  const bool& do_profile_wfs() const {
    return do_profile_wfs_;
  }

  void do_profile_wfs(const bool& do_profile_wfs_in) {
    do_profile_wfs_ = do_profile_wfs_in;
  }

  

  const bool& do_surface_wfs() const {
    return do_surface_wfs_;
  }

  void do_surface_wfs(const bool& do_surface_wfs_in) {
    do_surface_wfs_ = do_surface_wfs_in;
  }

  

  const bool& do_sleave_wfs() const {
    return do_sleave_wfs_;
  }

  void do_sleave_wfs(const bool& do_sleave_wfs_in) {
    do_sleave_wfs_ = do_sleave_wfs_in;
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

  

  const int& n_surface_wfs() const {
    return n_surface_wfs_;
  }

  void n_surface_wfs(const int& n_surface_wfs_in) {
    n_surface_wfs_ = n_surface_wfs_in;
  }

  

  const int& n_sleave_wfs() const {
    return n_sleave_wfs_;
  }

  void n_sleave_wfs(const int& n_sleave_wfs_in) {
    n_sleave_wfs_ = n_sleave_wfs_in;
  }

  

  const blitz::Array<double, 2>& lssl_slterm_isotropic() const {
    return lssl_slterm_isotropic_;
  }

  void lssl_slterm_isotropic(const blitz::Array<double, 2>& lssl_slterm_isotropic_in) {
    lssl_slterm_isotropic_ = lssl_slterm_isotropic_in;
  }

  

  const blitz::Array<double, 3>& lssl_slterm_f_0() const {
    return lssl_slterm_f_0_;
  }

  void lssl_slterm_f_0(const blitz::Array<double, 3>& lssl_slterm_f_0_in) {
    lssl_slterm_f_0_ = lssl_slterm_f_0_in;
  }

  

  const blitz::Array<double, 2>& l_deltau_input() const {
    return l_deltau_input_;
  }

  void l_deltau_input(const blitz::Array<double, 2>& l_deltau_input_in) {
    l_deltau_input_ = l_deltau_input_in;
  }

  

  const blitz::Array<double, 2>& l_omega_input() const {
    return l_omega_input_;
  }

  void l_omega_input(const blitz::Array<double, 2>& l_omega_input_in) {
    l_omega_input_ = l_omega_input_in;
  }

  

  const blitz::Array<double, 2>& l_asymm_input() const {
    return l_asymm_input_;
  }

  void l_asymm_input(const blitz::Array<double, 2>& l_asymm_input_in) {
    l_asymm_input_ = l_asymm_input_in;
  }

  

  const blitz::Array<double, 2>& l_d2s_scaling() const {
    return l_d2s_scaling_;
  }

  void l_d2s_scaling(const blitz::Array<double, 2>& l_d2s_scaling_in) {
    l_d2s_scaling_ = l_d2s_scaling_in;
  }

  
  // MANUAL EDIT //
  blitz::Array<double, 3>& ls_brdf_f_0() {
    return ls_brdf_f_0_;
  }
  // MANUAL EDIT //

  const blitz::Array<double, 3>& ls_brdf_f_0() const {
    return ls_brdf_f_0_;
  }

  void ls_brdf_f_0(const blitz::Array<double, 3>& ls_brdf_f_0_in) {
    ls_brdf_f_0_ = ls_brdf_f_0_in;
  }

  // MANUAL EDIT //
  blitz::Array<double, 2>& ls_brdf_f() {
    return ls_brdf_f_;
  }
  // MANUAL EDIT //

  const blitz::Array<double, 2>& ls_brdf_f() const {
    return ls_brdf_f_;
  }

  void ls_brdf_f(const blitz::Array<double, 2>& ls_brdf_f_in) {
    ls_brdf_f_ = ls_brdf_f_in;
  }

  // MANUAL EDIT //
  blitz::Array<double, 3>& ls_ubrdf_f() {
    return ls_ubrdf_f_;
  }
  // MANUAL EDIT //

  const blitz::Array<double, 3>& ls_ubrdf_f() const {
    return ls_ubrdf_f_;
  }

  void ls_ubrdf_f(const blitz::Array<double, 3>& ls_ubrdf_f_in) {
    ls_ubrdf_f_ = ls_ubrdf_f_in;
  }

  
  // MANUAL EDIT //
  blitz::Array<double, 1>& ls_emissivity() {
    return ls_emissivity_;
  }
  // MANUAL EDIT //

  const blitz::Array<double, 1>& ls_emissivity() const {
    return ls_emissivity_;
  }

  void ls_emissivity(const blitz::Array<double, 1>& ls_emissivity_in) {
    ls_emissivity_ = ls_emissivity_in;
  }

  

  const blitz::Array<double, 1>& intensity_toa() const {
    return intensity_toa_;
  }

  void intensity_toa(const blitz::Array<double, 1>& intensity_toa_in) {
    intensity_toa_ = intensity_toa_in;
  }

  

  const blitz::Array<double, 3>& profilewf_toa() const {
    return profilewf_toa_;
  }

  void profilewf_toa(const blitz::Array<double, 3>& profilewf_toa_in) {
    profilewf_toa_ = profilewf_toa_in;
  }

  

  const blitz::Array<double, 2>& surfacewf_toa() const {
    return surfacewf_toa_;
  }

  void surfacewf_toa(const blitz::Array<double, 2>& surfacewf_toa_in) {
    surfacewf_toa_ = surfacewf_toa_in;
  }

  

  const blitz::Array<double, 1>& intensity_boa() const {
    return intensity_boa_;
  }

  void intensity_boa(const blitz::Array<double, 1>& intensity_boa_in) {
    intensity_boa_ = intensity_boa_in;
  }

  

  const blitz::Array<double, 3>& profilewf_boa() const {
    return profilewf_boa_;
  }

  void profilewf_boa(const blitz::Array<double, 3>& profilewf_boa_in) {
    profilewf_boa_ = profilewf_boa_in;
  }

  

  const blitz::Array<double, 2>& surfacewf_boa() const {
    return surfacewf_boa_;
  }

  void surfacewf_boa(const blitz::Array<double, 2>& surfacewf_boa_in) {
    surfacewf_boa_ = surfacewf_boa_in;
  }

  

  const blitz::Array<double, 2>& radlevel_up() const {
    return radlevel_up_;
  }

  void radlevel_up(const blitz::Array<double, 2>& radlevel_up_in) {
    radlevel_up_ = radlevel_up_in;
  }

  

  const blitz::Array<double, 2>& radlevel_dn() const {
    return radlevel_dn_;
  }

  void radlevel_dn(const blitz::Array<double, 2>& radlevel_dn_in) {
    radlevel_dn_ = radlevel_dn_in;
  }

  

  const int& n_geometries() const {
    return n_geometries_;
  }

  

  const blitz::Array<double, 4>& profjaclevel_up() const {
    return profjaclevel_up_;
  }

  void profjaclevel_up(const blitz::Array<double, 4>& profjaclevel_up_in) {
    profjaclevel_up_ = profjaclevel_up_in;
  }

  

  const blitz::Array<double, 4>& profjaclevel_dn() const {
    return profjaclevel_dn_;
  }

  void profjaclevel_dn(const blitz::Array<double, 4>& profjaclevel_dn_in) {
    profjaclevel_dn_ = profjaclevel_dn_in;
  }

  

  const blitz::Array<double, 3>& surfjaclevel_up() const {
    return surfjaclevel_up_;
  }

  void surfjaclevel_up(const blitz::Array<double, 3>& surfjaclevel_up_in) {
    surfjaclevel_up_ = surfjaclevel_up_in;
  }

  

  const blitz::Array<double, 3>& surfjaclevel_dn() const {
    return surfjaclevel_dn_;
  }

  void surfjaclevel_dn(const blitz::Array<double, 3>& surfjaclevel_dn_in) {
    surfjaclevel_dn_ = surfjaclevel_dn_in;
  }

  

  const blitz::Array<double, 2>& fluxes_toa() const {
    return fluxes_toa_;
  }

  void fluxes_toa(const blitz::Array<double, 2>& fluxes_toa_in) {
    fluxes_toa_ = fluxes_toa_in;
  }

  

  const blitz::Array<double, 4>& profjacfluxes_toa() const {
    return profjacfluxes_toa_;
  }

  void profjacfluxes_toa(const blitz::Array<double, 4>& profjacfluxes_toa_in) {
    profjacfluxes_toa_ = profjacfluxes_toa_in;
  }

  

  const blitz::Array<double, 3>& surfjacfluxes_toa() const {
    return surfjacfluxes_toa_;
  }

  void surfjacfluxes_toa(const blitz::Array<double, 3>& surfjacfluxes_toa_in) {
    surfjacfluxes_toa_ = surfjacfluxes_toa_in;
  }

  

  const blitz::Array<double, 2>& fluxes_boa() const {
    return fluxes_boa_;
  }

  void fluxes_boa(const blitz::Array<double, 2>& fluxes_boa_in) {
    fluxes_boa_ = fluxes_boa_in;
  }

  

  const blitz::Array<double, 4>& profjacfluxes_boa() const {
    return profjacfluxes_boa_;
  }

  void profjacfluxes_boa(const blitz::Array<double, 4>& profjacfluxes_boa_in) {
    profjacfluxes_boa_ = profjacfluxes_boa_in;
  }

  

  const blitz::Array<double, 3>& surfjacfluxes_boa() const {
    return surfjacfluxes_boa_;
  }

  void surfjacfluxes_boa(const blitz::Array<double, 3>& surfjacfluxes_boa_in) {
    surfjacfluxes_boa_ = surfjacfluxes_boa_in;
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
    
    twostream_lps_master_m_twostream_lps_master_wrap(&maxlayers_, &maxtotal_, &maxmessages_, &maxbeams_, &max_geometries_, &max_user_streams_, &max_user_relazms_, &max_user_obsgeoms_, &max_atmoswfs_, &max_surfacewfs_, &max_sleavewfs_, &do_upwelling_, &do_dnwelling_, &do_plane_parallel_, &do_2s_levelout_, &do_mvout_only_, &do_additional_mvout_, &do_solar_sources_, &do_thermal_emission_, &do_surface_emission_, &do_d2s_scaling_, &do_brdf_surface_, &do_user_obsgeoms_, &do_surface_leaving_, &do_sl_isotropic_, &do_pentadiag_inverse_, &bvpindex_, &bvpscalefactor_, &taylor_order_, &taylor_small_, &tcutoff_, &nlayers_, &ntotal_, &stream_value_, &n_user_obsgeoms_, user_obsgeoms_.dataFirst(), &n_user_streams_, user_angles_.dataFirst(), &n_user_relazms_, user_relazms_.dataFirst(), &flux_factor_, &nbeams_, beam_szas_.dataFirst(), &earth_radius_, height_grid_.dataFirst(), deltau_input_.dataFirst(), omega_input_.dataFirst(), asymm_input_.dataFirst(), d2s_scaling_.dataFirst(), thermal_bb_input_.dataFirst(), &lambertian_albedo_, brdf_f_0_.dataFirst(), brdf_f_.dataFirst(), ubrdf_f_.dataFirst(), &emissivity_, &surfbb_, slterm_isotropic_.dataFirst(), slterm_f_0_.dataFirst(), &do_profile_wfs_, &do_surface_wfs_, &do_sleave_wfs_, layer_vary_flag_.dataFirst(), layer_vary_number_.dataFirst(), &n_surface_wfs_, &n_sleave_wfs_, lssl_slterm_isotropic_.dataFirst(), lssl_slterm_f_0_.dataFirst(), l_deltau_input_.dataFirst(), l_omega_input_.dataFirst(), l_asymm_input_.dataFirst(), l_d2s_scaling_.dataFirst(), ls_brdf_f_0_.dataFirst(), ls_brdf_f_.dataFirst(), ls_ubrdf_f_.dataFirst(), ls_emissivity_.dataFirst(), intensity_toa_.dataFirst(), profilewf_toa_.dataFirst(), surfacewf_toa_.dataFirst(), intensity_boa_.dataFirst(), profilewf_boa_.dataFirst(), surfacewf_boa_.dataFirst(), radlevel_up_.dataFirst(), radlevel_dn_.dataFirst(), &n_geometries_, profjaclevel_up_.dataFirst(), profjaclevel_dn_.dataFirst(), surfjaclevel_up_.dataFirst(), surfjaclevel_dn_.dataFirst(), fluxes_toa_.dataFirst(), profjacfluxes_toa_.dataFirst(), surfjacfluxes_toa_.dataFirst(), fluxes_boa_.dataFirst(), profjacfluxes_boa_.dataFirst(), surfjacfluxes_boa_.dataFirst(), &status_inputcheck_, &c_nmessages_, &c_messages_shape_1, &c_messages_len, c_messages_.dataFirst(), &c_actions_shape_1, &c_actions_len, c_actions_.dataFirst(), &status_execution_, &e_message_len, e_message_.dataFirst(), &e_trace_1_len, e_trace_1_.dataFirst(), &e_trace_2_len, e_trace_2_.dataFirst());
    
  }

  // MANUAL CHANGE
  void print(std::ostream &output_stream) const {
    output_stream << "Twostream_Lps_Master:" << std::endl
      << "            maxlayers: " << maxlayers()  << std::endl
      << "             maxtotal: " << maxtotal()  << std::endl
      << "          maxmessages: " << maxmessages()  << std::endl
      << "             maxbeams: " << maxbeams()  << std::endl
      << "       max_geometries: " << max_geometries()  << std::endl
      << "     max_user_streams: " << max_user_streams()  << std::endl
      << "     max_user_relazms: " << max_user_relazms()  << std::endl
      << "    max_user_obsgeoms: " << max_user_obsgeoms()  << std::endl
      << "         max_atmoswfs: " << max_atmoswfs()  << std::endl
      << "       max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "        max_sleavewfs: " << max_sleavewfs()  << std::endl
      << "         do_upwelling: " << do_upwelling()  << std::endl
      << "         do_dnwelling: " << do_dnwelling()  << std::endl
      << "    do_plane_parallel: " << do_plane_parallel()  << std::endl
      << "       do_2s_levelout: " << do_2s_levelout()  << std::endl
      << "        do_mvout_only: " << do_mvout_only()  << std::endl
      << "  do_additional_mvout: " << do_additional_mvout()  << std::endl
      << "     do_solar_sources: " << do_solar_sources()  << std::endl
      << "  do_thermal_emission: " << do_thermal_emission()  << std::endl
      << "  do_surface_emission: " << do_surface_emission()  << std::endl
      << "       do_d2s_scaling: " << do_d2s_scaling()  << std::endl
      << "      do_brdf_surface: " << do_brdf_surface()  << std::endl
      << "     do_user_obsgeoms: " << do_user_obsgeoms()  << std::endl
      << "   do_surface_leaving: " << do_surface_leaving()  << std::endl
      << "      do_sl_isotropic: " << do_sl_isotropic()  << std::endl
      << " do_pentadiag_inverse: " << do_pentadiag_inverse()  << std::endl
      << "             bvpindex: " << bvpindex()  << std::endl
      << "       bvpscalefactor: " << bvpscalefactor()  << std::endl
      << "         taylor_order: " << taylor_order()  << std::endl
      << "         taylor_small: " << taylor_small()  << std::endl
      << "              tcutoff: " << tcutoff()  << std::endl
      << "              nlayers: " << nlayers()  << std::endl
      << "               ntotal: " << ntotal()  << std::endl
      << "         stream_value: " << stream_value()  << std::endl
      << "      n_user_obsgeoms: " << n_user_obsgeoms()  << std::endl
      << "        user_obsgeoms: " << std::endl << user_obsgeoms()  << std::endl
      << "       n_user_streams: " << n_user_streams()  << std::endl
      << "          user_angles: " << std::endl << user_angles()  << std::endl
      << "       n_user_relazms: " << n_user_relazms()  << std::endl
      << "         user_relazms: " << std::endl << user_relazms()  << std::endl
      << "          flux_factor: " << flux_factor()  << std::endl
      << "               nbeams: " << nbeams()  << std::endl
      << "            beam_szas: " << std::endl << beam_szas()  << std::endl
      << "         earth_radius: " << earth_radius()  << std::endl
      << "          height_grid: " << std::endl << height_grid()  << std::endl
      << "         deltau_input: " << std::endl << deltau_input()  << std::endl
      << "          omega_input: " << std::endl << omega_input()  << std::endl
      << "          asymm_input: " << std::endl << asymm_input()  << std::endl
      << "          d2s_scaling: " << std::endl << d2s_scaling()  << std::endl
      << "     thermal_bb_input: " << std::endl << thermal_bb_input()  << std::endl
      << "    lambertian_albedo: " << lambertian_albedo()  << std::endl
      << "             brdf_f_0: " << std::endl << brdf_f_0()  << std::endl
      << "               brdf_f: " << std::endl << brdf_f()  << std::endl
      << "              ubrdf_f: " << std::endl << ubrdf_f()  << std::endl
      << "           emissivity: " << emissivity()  << std::endl
      << "               surfbb: " << surfbb()  << std::endl
      << "     slterm_isotropic: " << std::endl << slterm_isotropic()  << std::endl
      << "           slterm_f_0: " << std::endl << slterm_f_0()  << std::endl
      << "       do_profile_wfs: " << do_profile_wfs()  << std::endl
      << "       do_surface_wfs: " << do_surface_wfs()  << std::endl
      << "        do_sleave_wfs: " << do_sleave_wfs()  << std::endl
      << "      layer_vary_flag: " << std::endl << layer_vary_flag()  << std::endl
      << "    layer_vary_number: " << std::endl << layer_vary_number()  << std::endl
      << "        n_surface_wfs: " << n_surface_wfs()  << std::endl
      << "         n_sleave_wfs: " << n_sleave_wfs()  << std::endl
      << "lssl_slterm_isotropic: " << std::endl << lssl_slterm_isotropic()  << std::endl
      << "      lssl_slterm_f_0: " << std::endl << lssl_slterm_f_0()  << std::endl
      << "       l_deltau_input: " << std::endl << l_deltau_input()  << std::endl
      << "        l_omega_input: " << std::endl << l_omega_input()  << std::endl
      << "        l_asymm_input: " << std::endl << l_asymm_input()  << std::endl
      << "        l_d2s_scaling: " << std::endl << l_d2s_scaling()  << std::endl
      << "          ls_brdf_f_0: " << std::endl << ls_brdf_f_0()  << std::endl
      << "            ls_brdf_f: " << std::endl << ls_brdf_f()  << std::endl
      << "           ls_ubrdf_f: " << std::endl << ls_ubrdf_f()  << std::endl
      << "        ls_emissivity: " << std::endl << ls_emissivity()  << std::endl
      << "        intensity_toa: " << std::endl << intensity_toa()  << std::endl
      << "        profilewf_toa: " << std::endl << profilewf_toa()  << std::endl
      << "        surfacewf_toa: " << std::endl << surfacewf_toa()  << std::endl
      << "        intensity_boa: " << std::endl << intensity_boa()  << std::endl
      << "        profilewf_boa: " << std::endl << profilewf_boa()  << std::endl
      << "        surfacewf_boa: " << std::endl << surfacewf_boa()  << std::endl
      << "          radlevel_up: " << std::endl << radlevel_up()  << std::endl
      << "          radlevel_dn: " << std::endl << radlevel_dn()  << std::endl
      << "         n_geometries: " << n_geometries()  << std::endl
      << "      profjaclevel_up: " << std::endl << profjaclevel_up()  << std::endl
      << "      profjaclevel_dn: " << std::endl << profjaclevel_dn()  << std::endl
      << "      surfjaclevel_up: " << std::endl << surfjaclevel_up()  << std::endl
      << "      surfjaclevel_dn: " << std::endl << surfjaclevel_dn()  << std::endl
      << "           fluxes_toa: " << std::endl << fluxes_toa()  << std::endl
      << "    profjacfluxes_toa: " << std::endl << profjacfluxes_toa()  << std::endl
      << "    surfjacfluxes_toa: " << std::endl << surfjacfluxes_toa()  << std::endl
      << "           fluxes_boa: " << std::endl << fluxes_boa()  << std::endl
      << "    profjacfluxes_boa: " << std::endl << profjacfluxes_boa()  << std::endl
      << "    surfjacfluxes_boa: " << std::endl << surfjacfluxes_boa()  << std::endl
      << "    status_inputcheck: " << status_inputcheck()  << std::endl
      << "          c_nmessages: " << c_nmessages()  << std::endl
      << "           c_messages: " << std::endl;
    std::vector< std::string > c_messages_lcl = c_messages();
    for(unsigned int idx = 0; idx < c_messages_lcl.size(); idx++)
      if ( c_messages_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << c_messages_lcl[idx] << "\"" << std::endl;
    output_stream
      << "            c_actions: " << std::endl;
    std::vector< std::string > c_actions_lcl = c_actions();
    for(unsigned int idx = 0; idx < c_actions_lcl.size(); idx++)
      if ( c_actions_lcl[idx].length() > 0 )
        output_stream << "  [" << idx << "]: \"" << c_actions_lcl[idx] << "\"" << std::endl;
    output_stream
      << "     status_execution: " << status_execution()  << std::endl
      << "            e_message: " << "\"" << e_message() << "\"" << std::endl
      << "            e_trace_1: " << "\"" << e_trace_1() << "\"" << std::endl
      << "            e_trace_2: " << "\"" << e_trace_2() << "\"" << std::endl;
  }
  // MANUAL CHANGE

private:
  int maxlayers_;
  int maxtotal_;
  int maxmessages_;
  int maxbeams_;
  int max_geometries_;
  int max_user_streams_;
  int max_user_relazms_;
  int max_user_obsgeoms_;
  int max_atmoswfs_;
  int max_surfacewfs_;
  int max_sleavewfs_;
  bool do_upwelling_;
  bool do_dnwelling_;
  bool do_plane_parallel_;
  bool do_2s_levelout_;
  bool do_mvout_only_;
  bool do_additional_mvout_;
  bool do_solar_sources_;
  bool do_thermal_emission_;
  bool do_surface_emission_;
  bool do_d2s_scaling_;
  bool do_brdf_surface_;
  bool do_user_obsgeoms_;
  bool do_surface_leaving_;
  bool do_sl_isotropic_;
  bool do_pentadiag_inverse_;
  int bvpindex_;
  double bvpscalefactor_;
  int taylor_order_;
  double taylor_small_;
  double tcutoff_;
  int nlayers_;
  int ntotal_;
  double stream_value_;
  int n_user_obsgeoms_;
  blitz::Array<double, 2> user_obsgeoms_;
  int n_user_streams_;
  blitz::Array<double, 1> user_angles_;
  int n_user_relazms_;
  blitz::Array<double, 1> user_relazms_;
  double flux_factor_;
  int nbeams_;
  blitz::Array<double, 1> beam_szas_;
  double earth_radius_;
  blitz::Array<double, 1> height_grid_;
  blitz::Array<double, 1> deltau_input_;
  blitz::Array<double, 1> omega_input_;
  blitz::Array<double, 1> asymm_input_;
  blitz::Array<double, 1> d2s_scaling_;
  blitz::Array<double, 1> thermal_bb_input_;
  double lambertian_albedo_;
  blitz::Array<double, 2> brdf_f_0_;
  blitz::Array<double, 1> brdf_f_;
  blitz::Array<double, 2> ubrdf_f_;
  double emissivity_;
  double surfbb_;
  blitz::Array<double, 1> slterm_isotropic_;
  blitz::Array<double, 2> slterm_f_0_;
  bool do_profile_wfs_;
  bool do_surface_wfs_;
  bool do_sleave_wfs_;
  blitz::Array<bool, 1> layer_vary_flag_;
  blitz::Array<int, 1> layer_vary_number_;
  int n_surface_wfs_;
  int n_sleave_wfs_;
  blitz::Array<double, 2> lssl_slterm_isotropic_;
  blitz::Array<double, 3> lssl_slterm_f_0_;
  blitz::Array<double, 2> l_deltau_input_;
  blitz::Array<double, 2> l_omega_input_;
  blitz::Array<double, 2> l_asymm_input_;
  blitz::Array<double, 2> l_d2s_scaling_;
  blitz::Array<double, 3> ls_brdf_f_0_;
  blitz::Array<double, 2> ls_brdf_f_;
  blitz::Array<double, 3> ls_ubrdf_f_;
  blitz::Array<double, 1> ls_emissivity_;
  blitz::Array<double, 1> intensity_toa_;
  blitz::Array<double, 3> profilewf_toa_;
  blitz::Array<double, 2> surfacewf_toa_;
  blitz::Array<double, 1> intensity_boa_;
  blitz::Array<double, 3> profilewf_boa_;
  blitz::Array<double, 2> surfacewf_boa_;
  blitz::Array<double, 2> radlevel_up_;
  blitz::Array<double, 2> radlevel_dn_;
  int n_geometries_;
  blitz::Array<double, 4> profjaclevel_up_;
  blitz::Array<double, 4> profjaclevel_dn_;
  blitz::Array<double, 3> surfjaclevel_up_;
  blitz::Array<double, 3> surfjaclevel_dn_;
  blitz::Array<double, 2> fluxes_toa_;
  blitz::Array<double, 4> profjacfluxes_toa_;
  blitz::Array<double, 3> surfjacfluxes_toa_;
  blitz::Array<double, 2> fluxes_boa_;
  blitz::Array<double, 4> profjacfluxes_boa_;
  blitz::Array<double, 3> surfjacfluxes_boa_;
  int status_inputcheck_;
  int c_nmessages_;
  blitz::Array<char, 2> c_messages_;
  blitz::Array<char, 2> c_actions_;
  int status_execution_;
  blitz::Array<char, 1> e_message_;
  blitz::Array<char, 1> e_trace_1_;
  blitz::Array<char, 1> e_trace_2_;
  // MANUAL CHANGE
  Twostream_Lps_Master() {}
  // MANUAL CHANGE
};

}

#endif
