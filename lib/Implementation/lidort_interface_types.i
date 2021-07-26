// This file was auto-generated

%include "common.i"

%{
#include "lidort_interface_types.h"
%}

%fp_shared_ptr(FullPhysics::Lidort_Pars);

namespace FullPhysics {

%nodefaultctor Lidort_Pars;
%nodefaultctor Lidort_Structure;


struct Lidort_Pars {

  const char lidort_version_number[5];
  const int lidort_inunit;
  const int lidort_scenunit;
  const int lidort_funit;
  const int lidort_resunit;
  const int lidort_errunit;
  const int lidort_dbgunit;
  const int max_messages;
  const int maxstreams;
  const int maxlayers;
  const int maxfinelayers;
  const int maxmoments_input;
  const int max_thermal_coeffs;
  const int maxbeams;
  const int max_user_streams;
  const int max_user_relazms;
  const int max_user_obsgeoms;
  const int max_user_levels;
  const int max_partlayers;
  const int max_taylor_terms;
  const int max_directions;
  const int max_brdf_kernels;
  const int max_brdf_parameters;
  const int maxstreams_brdf;
  const int max_msrs_muquad;
  const int max_msrs_phiquad;
  const int maxstreams_scaling;
  const int max_atmoswfs;
  const int max_surfacewfs;
  const int max_sleavewfs;
  const int max_geometries;
  const int max_allstrms;
  const int max_allstrms_p1;
  const int maxmoments;
  const int maxfourier;
  const int maxsthalf_brdf;
  const int maxstreams_2;
  const int maxstreams_p1;
  const int maxtotal;
  const int maxbandtotal;
  const double one;
  const double zero;
  const double onep5;
  const double two;
  const double three;
  const double four;
  const double quarter;
  const double half;
  const double minus_one;
  const double minus_two;
  const double pie;
  const double deg_to_rad;
  const double pi2;
  const double pi4;
  const double pio2;
  const double pio4;
  const double eps3;
  const double eps4;
  const double eps5;
  const double taylor_small;
  const double smallnum;
  const double bigexp;
  const double hopital_tolerance;
  const double omega_smallnum;
  const double max_tau_spath;
  const double max_tau_upath;
  const double max_tau_qpath;
  const int lidort_serious;
  const int lidort_warning;
  const int lidort_info;
  const int lidort_debug;
  const int lidort_success;
  const int upidx;
  const int dnidx;
  const int lambertian_idx;
  const int rossthin_idx;
  const int rossthick_idx;
  const int lisparse_idx;
  const int lidense_idx;
  const int hapke_idx;
  const int roujean_idx;
  const int rahman_idx;
  const int coxmunk_idx;
  const int bpdfsoil_idx;
  const int bpdfvegn_idx;
  const int bpdfndvi_idx;
  const int newcmglint_idx;
  const int rtkhotspot_idx;
  const int modfresnel_idx;
  const int snowbrdf_idx;
  const int maxbrdf_idx;
  
  static Lidort_Pars& instance();

};


class Lidort_Structure {
public:
  void* fortran_type_ptr();

  std::string print_to_string() const;
};

class Brdf_Linsup_Inputs : public Lidort_Structure {
public:
  Brdf_Linsup_Inputs();
  Brdf_Linsup_Inputs(const Brdf_Linsup_Inputs& src);
  ~Brdf_Linsup_Inputs();

  const blitz::Array<bool, 1> bs_do_kernel_factor_wfs() const;
  void bs_do_kernel_factor_wfs(const blitz::Array<bool, 1>& bs_do_kernel_factor_wfs_in);
  
  const blitz::Array<bool, 2> bs_do_kernel_params_wfs() const;
  void bs_do_kernel_params_wfs(const blitz::Array<bool, 2>& bs_do_kernel_params_wfs_in);
  
  const blitz::Array<bool, 1> bs_do_kparams_derivs() const;
  void bs_do_kparams_derivs(const blitz::Array<bool, 1>& bs_do_kparams_derivs_in);
  
  const int& bs_n_surface_wfs() const;
  void bs_n_surface_wfs(const int& bs_n_surface_wfs_in);
  
  const int& bs_n_kernel_factor_wfs() const;
  void bs_n_kernel_factor_wfs(const int& bs_n_kernel_factor_wfs_in);
  
  const int& bs_n_kernel_params_wfs() const;
  void bs_n_kernel_params_wfs(const int& bs_n_kernel_params_wfs_in);
  
  const bool bs_do_bsavalue_wf() const;
  void bs_do_bsavalue_wf(const bool& bs_do_bsavalue_wf_in);
  
  const bool bs_do_wsavalue_wf() const;
  void bs_do_wsavalue_wf(const bool& bs_do_wsavalue_wf_in);
  
  const bool bs_do_windspeed_wf() const;
  void bs_do_windspeed_wf(const bool& bs_do_windspeed_wf_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Brdf_Linsup_Outputs : public Lidort_Structure {
public:
  Brdf_Linsup_Outputs();
  Brdf_Linsup_Outputs(const Brdf_Linsup_Outputs& src);
  ~Brdf_Linsup_Outputs();

  const blitz::Array<double, 4>& bs_ls_dbounce_brdfunc() const;
  void bs_ls_dbounce_brdfunc(const blitz::Array<double, 4>& bs_ls_dbounce_brdfunc_in);
  
  const blitz::Array<double, 4>& bs_ls_brdf_f_0() const;
  void bs_ls_brdf_f_0(const blitz::Array<double, 4>& bs_ls_brdf_f_0_in);
  
  const blitz::Array<double, 4>& bs_ls_brdf_f() const;
  void bs_ls_brdf_f(const blitz::Array<double, 4>& bs_ls_brdf_f_in);
  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f_0() const;
  void bs_ls_user_brdf_f_0(const blitz::Array<double, 4>& bs_ls_user_brdf_f_0_in);
  
  const blitz::Array<double, 4>& bs_ls_user_brdf_f() const;
  void bs_ls_user_brdf_f(const blitz::Array<double, 4>& bs_ls_user_brdf_f_in);
  
  const blitz::Array<double, 2>& bs_ls_emissivity() const;
  void bs_ls_emissivity(const blitz::Array<double, 2>& bs_ls_emissivity_in);
  
  const blitz::Array<double, 2>& bs_ls_user_emissivity() const;
  void bs_ls_user_emissivity(const blitz::Array<double, 2>& bs_ls_user_emissivity_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Brdf_Sup_Inputs : public Lidort_Structure {
public:
  Brdf_Sup_Inputs();
  Brdf_Sup_Inputs(const Brdf_Sup_Inputs& src);
  ~Brdf_Sup_Inputs();

  const bool bs_do_brdf_surface() const;
  void bs_do_brdf_surface(const bool& bs_do_brdf_surface_in);
  
  const bool bs_do_surface_emission() const;
  void bs_do_surface_emission(const bool& bs_do_surface_emission_in);
  
  const bool bs_do_solar_sources() const;
  void bs_do_solar_sources(const bool& bs_do_solar_sources_in);
  
  const bool bs_do_user_streams() const;
  void bs_do_user_streams(const bool& bs_do_user_streams_in);
  
  const bool bs_do_user_obsgeoms() const;
  void bs_do_user_obsgeoms(const bool& bs_do_user_obsgeoms_in);
  
  const bool bs_do_doublet_geometry() const;
  void bs_do_doublet_geometry(const bool& bs_do_doublet_geometry_in);
  
  const int& bs_nstreams() const;
  void bs_nstreams(const int& bs_nstreams_in);
  
  const int& bs_nbeams() const;
  void bs_nbeams(const int& bs_nbeams_in);
  
  const blitz::Array<double, 1>& bs_beam_szas() const;
  void bs_beam_szas(const blitz::Array<double, 1>& bs_beam_szas_in);
  
  const int& bs_n_user_relazms() const;
  void bs_n_user_relazms(const int& bs_n_user_relazms_in);
  
  const blitz::Array<double, 1>& bs_user_relazms() const;
  void bs_user_relazms(const blitz::Array<double, 1>& bs_user_relazms_in);
  
  const int& bs_n_user_streams() const;
  void bs_n_user_streams(const int& bs_n_user_streams_in);
  
  const blitz::Array<double, 1>& bs_user_angles_input() const;
  void bs_user_angles_input(const blitz::Array<double, 1>& bs_user_angles_input_in);
  
  const int& bs_n_user_obsgeoms() const;
  void bs_n_user_obsgeoms(const int& bs_n_user_obsgeoms_in);
  
  const blitz::Array<double, 2>& bs_user_obsgeoms() const;
  void bs_user_obsgeoms(const blitz::Array<double, 2>& bs_user_obsgeoms_in);
  
  const int& bs_n_user_doublets() const;
  void bs_n_user_doublets(const int& bs_n_user_doublets_in);
  
  const blitz::Array<double, 2>& bs_user_doublets() const;
  void bs_user_doublets(const blitz::Array<double, 2>& bs_user_doublets_in);
  
  const int& bs_n_brdf_kernels() const;
  void bs_n_brdf_kernels(const int& bs_n_brdf_kernels_in);
  
  const std::vector< std::string > bs_brdf_names() const;
  
  const blitz::Array<int, 1>& bs_which_brdf() const;
  void bs_which_brdf(const blitz::Array<int, 1>& bs_which_brdf_in);
  
  const blitz::Array<int, 1>& bs_n_brdf_parameters() const;
  void bs_n_brdf_parameters(const blitz::Array<int, 1>& bs_n_brdf_parameters_in);
  
  const blitz::Array<double, 2>& bs_brdf_parameters() const;
  void bs_brdf_parameters(const blitz::Array<double, 2>& bs_brdf_parameters_in);
  
  const blitz::Array<bool, 1> bs_lambertian_kernel_flag() const;
  void bs_lambertian_kernel_flag(const blitz::Array<bool, 1>& bs_lambertian_kernel_flag_in);
  
  const blitz::Array<double, 1>& bs_brdf_factors() const;
  void bs_brdf_factors(const blitz::Array<double, 1>& bs_brdf_factors_in);
  
  const int& bs_nstreams_brdf() const;
  void bs_nstreams_brdf(const int& bs_nstreams_brdf_in);
  
  const bool bs_do_shadow_effect() const;
  void bs_do_shadow_effect(const bool& bs_do_shadow_effect_in);
  
  const bool bs_do_directbounce_only() const;
  void bs_do_directbounce_only(const bool& bs_do_directbounce_only_in);
  
  const bool bs_do_wsabsa_output() const;
  void bs_do_wsabsa_output(const bool& bs_do_wsabsa_output_in);
  
  const bool bs_do_wsa_scaling() const;
  void bs_do_wsa_scaling(const bool& bs_do_wsa_scaling_in);
  
  const bool bs_do_bsa_scaling() const;
  void bs_do_bsa_scaling(const bool& bs_do_bsa_scaling_in);
  
  const double& bs_wsa_value() const;
  void bs_wsa_value(const double& bs_wsa_value_in);
  
  const double& bs_bsa_value() const;
  void bs_bsa_value(const double& bs_bsa_value_in);
  
  const bool bs_do_newcmglint() const;
  void bs_do_newcmglint(const bool& bs_do_newcmglint_in);
  
  const double& bs_salinity() const;
  void bs_salinity(const double& bs_salinity_in);
  
  const double& bs_wavelength() const;
  void bs_wavelength(const double& bs_wavelength_in);
  
  const double& bs_windspeed() const;
  void bs_windspeed(const double& bs_windspeed_in);
  
  const blitz::Array<double, 1>& bs_winddir() const;
  void bs_winddir(const blitz::Array<double, 1>& bs_winddir_in);
  
  const bool bs_do_glintshadow() const;
  void bs_do_glintshadow(const bool& bs_do_glintshadow_in);
  
  const bool bs_do_foamoption() const;
  void bs_do_foamoption(const bool& bs_do_foamoption_in);
  
  const bool bs_do_facetisotropy() const;
  void bs_do_facetisotropy(const bool& bs_do_facetisotropy_in);
  
  const bool bs_do_glitter_msrcorr() const;
  void bs_do_glitter_msrcorr(const bool& bs_do_glitter_msrcorr_in);
  
  const bool bs_do_glitter_msrcorr_dbonly() const;
  void bs_do_glitter_msrcorr_dbonly(const bool& bs_do_glitter_msrcorr_dbonly_in);
  
  const int& bs_glitter_msrcorr_order() const;
  void bs_glitter_msrcorr_order(const int& bs_glitter_msrcorr_order_in);
  
  const int& bs_glitter_msrcorr_nmuquad() const;
  void bs_glitter_msrcorr_nmuquad(const int& bs_glitter_msrcorr_nmuquad_in);
  
  const int& bs_glitter_msrcorr_nphiquad() const;
  void bs_glitter_msrcorr_nphiquad(const int& bs_glitter_msrcorr_nphiquad_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Brdf_Sup_Outputs : public Lidort_Structure {
public:
  Brdf_Sup_Outputs();
  Brdf_Sup_Outputs(const Brdf_Sup_Outputs& src);
  ~Brdf_Sup_Outputs();

  const blitz::Array<double, 3>& bs_dbounce_brdfunc() const;
  void bs_dbounce_brdfunc(const blitz::Array<double, 3>& bs_dbounce_brdfunc_in);
  
  const blitz::Array<double, 3>& bs_brdf_f_0() const;
  void bs_brdf_f_0(const blitz::Array<double, 3>& bs_brdf_f_0_in);
  
  const blitz::Array<double, 3>& bs_brdf_f() const;
  void bs_brdf_f(const blitz::Array<double, 3>& bs_brdf_f_in);
  
  const blitz::Array<double, 3>& bs_user_brdf_f_0() const;
  void bs_user_brdf_f_0(const blitz::Array<double, 3>& bs_user_brdf_f_0_in);
  
  const blitz::Array<double, 3>& bs_user_brdf_f() const;
  void bs_user_brdf_f(const blitz::Array<double, 3>& bs_user_brdf_f_in);
  
  const blitz::Array<double, 1>& bs_emissivity() const;
  void bs_emissivity(const blitz::Array<double, 1>& bs_emissivity_in);
  
  const blitz::Array<double, 1>& bs_user_emissivity() const;
  void bs_user_emissivity(const blitz::Array<double, 1>& bs_user_emissivity_in);
  
  const double& bs_wsa_calculated() const;
  void bs_wsa_calculated(const double& bs_wsa_calculated_in);
  
  const blitz::Array<double, 1>& bs_wsa_kernels() const;
  void bs_wsa_kernels(const blitz::Array<double, 1>& bs_wsa_kernels_in);
  
  const double& bs_bsa_calculated() const;
  void bs_bsa_calculated(const double& bs_bsa_calculated_in);
  
  const blitz::Array<double, 1>& bs_bsa_kernels() const;
  void bs_bsa_kernels(const blitz::Array<double, 1>& bs_bsa_kernels_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Brdf_Input_Exception_Handling : public Lidort_Structure {
public:
  Brdf_Input_Exception_Handling();
  Brdf_Input_Exception_Handling(const Brdf_Input_Exception_Handling& src);
  ~Brdf_Input_Exception_Handling();

  const int& bs_status_inputread() const;
  void bs_status_inputread(const int& bs_status_inputread_in);
  
  const int& bs_ninputmessages() const;
  void bs_ninputmessages(const int& bs_ninputmessages_in);
  
  const std::vector< std::string > bs_inputmessages() const;
  
  const std::vector< std::string > bs_inputactions() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Brdf_Output_Exception_Handling : public Lidort_Structure {
public:
  Brdf_Output_Exception_Handling();
  Brdf_Output_Exception_Handling(const Brdf_Output_Exception_Handling& src);
  ~Brdf_Output_Exception_Handling();

  const int& bs_status_output() const;
  void bs_status_output(const int& bs_status_output_in);
  
  const int& bs_noutputmessages() const;
  void bs_noutputmessages(const int& bs_noutputmessages_in);
  
  const std::vector< std::string > bs_outputmessages() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Sleave_Sup_Inputs : public Lidort_Structure {
public:
  Sleave_Sup_Inputs();
  Sleave_Sup_Inputs(const Sleave_Sup_Inputs& src);
  ~Sleave_Sup_Inputs();

  const bool sl_do_sleaving() const;
  void sl_do_sleaving(const bool& sl_do_sleaving_in);
  
  const bool sl_do_isotropic() const;
  void sl_do_isotropic(const bool& sl_do_isotropic_in);
  
  const bool sl_do_roughsurface() const;
  void sl_do_roughsurface(const bool& sl_do_roughsurface_in);
  
  const bool sl_do_exact() const;
  void sl_do_exact(const bool& sl_do_exact_in);
  
  const bool sl_do_exactonly() const;
  void sl_do_exactonly(const bool& sl_do_exactonly_in);
  
  const bool sl_do_fluorescence() const;
  void sl_do_fluorescence(const bool& sl_do_fluorescence_in);
  
  const bool sl_do_solar_sources() const;
  void sl_do_solar_sources(const bool& sl_do_solar_sources_in);
  
  const std::string sl_sleave_datapath() const;
  
  const bool sl_do_user_streams() const;
  void sl_do_user_streams(const bool& sl_do_user_streams_in);
  
  const bool sl_do_user_obsgeoms() const;
  void sl_do_user_obsgeoms(const bool& sl_do_user_obsgeoms_in);
  
  const bool sl_do_doublet_geometry() const;
  void sl_do_doublet_geometry(const bool& sl_do_doublet_geometry_in);
  
  const int& sl_nstreams() const;
  void sl_nstreams(const int& sl_nstreams_in);
  
  const int& sl_nbeams() const;
  void sl_nbeams(const int& sl_nbeams_in);
  
  const blitz::Array<double, 1>& sl_beam_szas() const;
  void sl_beam_szas(const blitz::Array<double, 1>& sl_beam_szas_in);
  
  const int& sl_n_user_relazms() const;
  void sl_n_user_relazms(const int& sl_n_user_relazms_in);
  
  const blitz::Array<double, 1>& sl_user_relazms() const;
  void sl_user_relazms(const blitz::Array<double, 1>& sl_user_relazms_in);
  
  const int& sl_n_user_streams() const;
  void sl_n_user_streams(const int& sl_n_user_streams_in);
  
  const blitz::Array<double, 1>& sl_user_angles_input() const;
  void sl_user_angles_input(const blitz::Array<double, 1>& sl_user_angles_input_in);
  
  const int& sl_n_user_obsgeoms() const;
  void sl_n_user_obsgeoms(const int& sl_n_user_obsgeoms_in);
  
  const blitz::Array<double, 2>& sl_user_obsgeoms() const;
  void sl_user_obsgeoms(const blitz::Array<double, 2>& sl_user_obsgeoms_in);
  
  const int& sl_n_user_doublets() const;
  void sl_n_user_doublets(const int& sl_n_user_doublets_in);
  
  const blitz::Array<double, 2>& sl_user_doublets() const;
  void sl_user_doublets(const blitz::Array<double, 2>& sl_user_doublets_in);
  
  const double& sl_salinity() const;
  void sl_salinity(const double& sl_salinity_in);
  
  const double& sl_chlorconc() const;
  void sl_chlorconc(const double& sl_chlorconc_in);
  
  const double& sl_wavelength() const;
  void sl_wavelength(const double& sl_wavelength_in);
  
  const bool sl_azimuthdep() const;
  void sl_azimuthdep(const bool& sl_azimuthdep_in);
  
  const bool sl_do_fourier_output() const;
  void sl_do_fourier_output(const bool& sl_do_fourier_output_in);
  
  const double& sl_windspeed() const;
  void sl_windspeed(const double& sl_windspeed_in);
  
  const blitz::Array<double, 1>& sl_winddir() const;
  void sl_winddir(const blitz::Array<double, 1>& sl_winddir_in);
  
  const bool sl_do_glintshadow() const;
  void sl_do_glintshadow(const bool& sl_do_glintshadow_in);
  
  const bool sl_do_foamoption() const;
  void sl_do_foamoption(const bool& sl_do_foamoption_in);
  
  const bool sl_do_facetisotropy() const;
  void sl_do_facetisotropy(const bool& sl_do_facetisotropy_in);
  
  const double& sl_fl_wavelength() const;
  void sl_fl_wavelength(const double& sl_fl_wavelength_in);
  
  const double& sl_fl_latitude() const;
  void sl_fl_latitude(const double& sl_fl_latitude_in);
  
  const double& sl_fl_longitude() const;
  void sl_fl_longitude(const double& sl_fl_longitude_in);
  
  const blitz::Array<int, 1>& sl_fl_epoch() const;
  void sl_fl_epoch(const blitz::Array<int, 1>& sl_fl_epoch_in);
  
  const double& sl_fl_amplitude755() const;
  void sl_fl_amplitude755(const double& sl_fl_amplitude755_in);
  
  const bool sl_fl_do_datagaussian() const;
  void sl_fl_do_datagaussian(const bool& sl_fl_do_datagaussian_in);
  
  const blitz::Array<double, 2>& sl_fl_inputgaussians() const;
  void sl_fl_inputgaussians(const blitz::Array<double, 2>& sl_fl_inputgaussians_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Lincontrol : public Lidort_Structure {
public:
  Lidort_Fixed_Lincontrol();
  Lidort_Fixed_Lincontrol(const Lidort_Fixed_Lincontrol& src);
  ~Lidort_Fixed_Lincontrol();

  const blitz::Array<bool, 1> ts_layer_vary_flag() const;
  void ts_layer_vary_flag(const blitz::Array<bool, 1>& ts_layer_vary_flag_in);
  
  const blitz::Array<int, 1>& ts_layer_vary_number() const;
  void ts_layer_vary_number(const blitz::Array<int, 1>& ts_layer_vary_number_in);
  
  const int& ts_n_totalcolumn_wfs() const;
  void ts_n_totalcolumn_wfs(const int& ts_n_totalcolumn_wfs_in);
  
  const int& ts_n_surface_wfs() const;
  void ts_n_surface_wfs(const int& ts_n_surface_wfs_in);
  
  const int& ts_n_sleave_wfs() const;
  void ts_n_sleave_wfs(const int& ts_n_sleave_wfs_in);
  
  const std::vector< std::string > ts_columnwf_names() const;
  
  const std::vector< std::string > ts_profilewf_names() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Linoptical : public Lidort_Structure {
public:
  Lidort_Fixed_Linoptical();
  Lidort_Fixed_Linoptical(const Lidort_Fixed_Linoptical& src);
  ~Lidort_Fixed_Linoptical();

  const blitz::Array<double, 2>& ts_l_deltau_vert_input() const;
  void ts_l_deltau_vert_input(const blitz::Array<double, 2>& ts_l_deltau_vert_input_in);
  
  const blitz::Array<double, 2>& ts_l_omega_total_input() const;
  void ts_l_omega_total_input(const blitz::Array<double, 2>& ts_l_omega_total_input_in);
  
  const blitz::Array<double, 3>& ts_l_phasmoms_total_input() const;
  void ts_l_phasmoms_total_input(const blitz::Array<double, 3>& ts_l_phasmoms_total_input_in);
  
  const blitz::Array<double, 3>& ts_l_phasfunc_input_up() const;
  void ts_l_phasfunc_input_up(const blitz::Array<double, 3>& ts_l_phasfunc_input_up_in);
  
  const blitz::Array<double, 3>& ts_l_phasfunc_input_dn() const;
  void ts_l_phasfunc_input_dn(const blitz::Array<double, 3>& ts_l_phasfunc_input_dn_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Lininputs : public Lidort_Structure {
public:
  Lidort_Fixed_Lininputs();
  Lidort_Fixed_Lininputs(const Lidort_Fixed_Lininputs& src);
  ~Lidort_Fixed_Lininputs();

  const Lidort_Fixed_Lincontrol& cont() const;
  void cont(Lidort_Fixed_Lincontrol& cont_in);
  
  const Lidort_Fixed_Linoptical& optical() const;
  void optical(Lidort_Fixed_Linoptical& optical_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Lincontrol : public Lidort_Structure {
public:
  Lidort_Modified_Lincontrol();
  Lidort_Modified_Lincontrol(const Lidort_Modified_Lincontrol& src);
  ~Lidort_Modified_Lincontrol();

  const bool ts_do_column_linearization() const;
  void ts_do_column_linearization(const bool& ts_do_column_linearization_in);
  
  const bool ts_do_profile_linearization() const;
  void ts_do_profile_linearization(const bool& ts_do_profile_linearization_in);
  
  const bool ts_do_atmos_linearization() const;
  void ts_do_atmos_linearization(const bool& ts_do_atmos_linearization_in);
  
  const bool ts_do_surface_linearization() const;
  void ts_do_surface_linearization(const bool& ts_do_surface_linearization_in);
  
  const bool ts_do_linearization() const;
  void ts_do_linearization(const bool& ts_do_linearization_in);
  
  const bool ts_do_simulation_only() const;
  void ts_do_simulation_only(const bool& ts_do_simulation_only_in);
  
  const bool ts_do_atmos_lbbf() const;
  void ts_do_atmos_lbbf(const bool& ts_do_atmos_lbbf_in);
  
  const bool ts_do_surface_lbbf() const;
  void ts_do_surface_lbbf(const bool& ts_do_surface_lbbf_in);
  
  const bool ts_do_sleave_wfs() const;
  void ts_do_sleave_wfs(const bool& ts_do_sleave_wfs_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Lininputs : public Lidort_Structure {
public:
  Lidort_Modified_Lininputs();
  Lidort_Modified_Lininputs(const Lidort_Modified_Lininputs& src);
  ~Lidort_Modified_Lininputs();

  const Lidort_Modified_Lincontrol& mcont() const;
  void mcont(Lidort_Modified_Lincontrol& mcont_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linatmos : public Lidort_Structure {
public:
  Lidort_Linatmos();
  Lidort_Linatmos(const Lidort_Linatmos& src);
  ~Lidort_Linatmos();

  const blitz::Array<double, 4>& ts_columnwf() const;
  void ts_columnwf(const blitz::Array<double, 4>& ts_columnwf_in);
  
  const blitz::Array<double, 4>& ts_meani_diffuse_colwf() const;
  void ts_meani_diffuse_colwf(const blitz::Array<double, 4>& ts_meani_diffuse_colwf_in);
  
  const blitz::Array<double, 4>& ts_flux_diffuse_colwf() const;
  void ts_flux_diffuse_colwf(const blitz::Array<double, 4>& ts_flux_diffuse_colwf_in);
  
  const blitz::Array<double, 3>& ts_dnmeani_direct_colwf() const;
  void ts_dnmeani_direct_colwf(const blitz::Array<double, 3>& ts_dnmeani_direct_colwf_in);
  
  const blitz::Array<double, 3>& ts_dnflux_direct_colwf() const;
  void ts_dnflux_direct_colwf(const blitz::Array<double, 3>& ts_dnflux_direct_colwf_in);
  
  const blitz::Array<double, 5>& ts_profilewf() const;
  void ts_profilewf(const blitz::Array<double, 5>& ts_profilewf_in);
  
  const blitz::Array<double, 5>& ts_meani_diffuse_profwf() const;
  void ts_meani_diffuse_profwf(const blitz::Array<double, 5>& ts_meani_diffuse_profwf_in);
  
  const blitz::Array<double, 5>& ts_flux_diffuse_profwf() const;
  void ts_flux_diffuse_profwf(const blitz::Array<double, 5>& ts_flux_diffuse_profwf_in);
  
  const blitz::Array<double, 4>& ts_dnmeani_direct_profwf() const;
  void ts_dnmeani_direct_profwf(const blitz::Array<double, 4>& ts_dnmeani_direct_profwf_in);
  
  const blitz::Array<double, 4>& ts_dnflux_direct_profwf() const;
  void ts_dnflux_direct_profwf(const blitz::Array<double, 4>& ts_dnflux_direct_profwf_in);
  
  const blitz::Array<double, 4>& ts_abbwfs_jacobians() const;
  void ts_abbwfs_jacobians(const blitz::Array<double, 4>& ts_abbwfs_jacobians_in);
  
  const blitz::Array<double, 4>& ts_abbwfs_fluxes() const;
  void ts_abbwfs_fluxes(const blitz::Array<double, 4>& ts_abbwfs_fluxes_in);
  
  const blitz::Array<double, 3>& ts_albmed_user_profwf() const;
  void ts_albmed_user_profwf(const blitz::Array<double, 3>& ts_albmed_user_profwf_in);
  
  const blitz::Array<double, 3>& ts_trnmed_user_profwf() const;
  void ts_trnmed_user_profwf(const blitz::Array<double, 3>& ts_trnmed_user_profwf_in);
  
  const blitz::Array<double, 3>& ts_albmed_fluxes_profwf() const;
  void ts_albmed_fluxes_profwf(const blitz::Array<double, 3>& ts_albmed_fluxes_profwf_in);
  
  const blitz::Array<double, 3>& ts_trnmed_fluxes_profwf() const;
  void ts_trnmed_fluxes_profwf(const blitz::Array<double, 3>& ts_trnmed_fluxes_profwf_in);
  
  const blitz::Array<double, 3>& ts_transbeam_profwf() const;
  void ts_transbeam_profwf(const blitz::Array<double, 3>& ts_transbeam_profwf_in);
  
  const blitz::Array<double, 2>& ts_albmed_user_colwf() const;
  void ts_albmed_user_colwf(const blitz::Array<double, 2>& ts_albmed_user_colwf_in);
  
  const blitz::Array<double, 2>& ts_trnmed_user_colwf() const;
  void ts_trnmed_user_colwf(const blitz::Array<double, 2>& ts_trnmed_user_colwf_in);
  
  const blitz::Array<double, 2>& ts_albmed_fluxes_colwf() const;
  void ts_albmed_fluxes_colwf(const blitz::Array<double, 2>& ts_albmed_fluxes_colwf_in);
  
  const blitz::Array<double, 2>& ts_trnmed_fluxes_colwf() const;
  void ts_trnmed_fluxes_colwf(const blitz::Array<double, 2>& ts_trnmed_fluxes_colwf_in);
  
  const blitz::Array<double, 2>& ts_transbeam_colwf() const;
  void ts_transbeam_colwf(const blitz::Array<double, 2>& ts_transbeam_colwf_in);
  
  const blitz::Array<double, 3>& ts_planetary_transterm_profwf() const;
  void ts_planetary_transterm_profwf(const blitz::Array<double, 3>& ts_planetary_transterm_profwf_in);
  
  const blitz::Array<double, 2>& ts_planetary_sbterm_profwf() const;
  void ts_planetary_sbterm_profwf(const blitz::Array<double, 2>& ts_planetary_sbterm_profwf_in);
  
  const blitz::Array<double, 2>& ts_planetary_transterm_colwf() const;
  void ts_planetary_transterm_colwf(const blitz::Array<double, 2>& ts_planetary_transterm_colwf_in);
  
  const blitz::Array<double, 1>& ts_planetary_sbterm_colwf() const;
  void ts_planetary_sbterm_colwf(const blitz::Array<double, 1>& ts_planetary_sbterm_colwf_in);
  
  const blitz::Array<double, 3>& ts_lc_lostrans() const;
  void ts_lc_lostrans(const blitz::Array<double, 3>& ts_lc_lostrans_in);
  
  const blitz::Array<double, 3>& ts_lc_layer_mssts() const;
  void ts_lc_layer_mssts(const blitz::Array<double, 3>& ts_lc_layer_mssts_in);
  
  const blitz::Array<double, 2>& ts_lc_surf_mssts() const;
  void ts_lc_surf_mssts(const blitz::Array<double, 2>& ts_lc_surf_mssts_in);
  
  const blitz::Array<double, 3>& ts_lp_lostrans() const;
  void ts_lp_lostrans(const blitz::Array<double, 3>& ts_lp_lostrans_in);
  
  const blitz::Array<double, 4>& ts_lp_layer_mssts() const;
  void ts_lp_layer_mssts(const blitz::Array<double, 4>& ts_lp_layer_mssts_in);
  
  const blitz::Array<double, 3>& ts_lp_surf_mssts() const;
  void ts_lp_surf_mssts(const blitz::Array<double, 3>& ts_lp_surf_mssts_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsurf : public Lidort_Structure {
public:
  Lidort_Linsurf();
  Lidort_Linsurf(const Lidort_Linsurf& src);
  ~Lidort_Linsurf();

  const blitz::Array<double, 4>& ts_surfacewf() const;
  void ts_surfacewf(const blitz::Array<double, 4>& ts_surfacewf_in);
  
  const blitz::Array<double, 4>& ts_meani_diffuse_surfwf() const;
  void ts_meani_diffuse_surfwf(const blitz::Array<double, 4>& ts_meani_diffuse_surfwf_in);
  
  const blitz::Array<double, 4>& ts_flux_diffuse_surfwf() const;
  void ts_flux_diffuse_surfwf(const blitz::Array<double, 4>& ts_flux_diffuse_surfwf_in);
  
  const blitz::Array<double, 3>& ts_sbbwfs_jacobians() const;
  void ts_sbbwfs_jacobians(const blitz::Array<double, 3>& ts_sbbwfs_jacobians_in);
  
  const blitz::Array<double, 3>& ts_sbbwfs_fluxes() const;
  void ts_sbbwfs_fluxes(const blitz::Array<double, 3>& ts_sbbwfs_fluxes_in);
  
  const blitz::Array<double, 3>& ts_ls_layer_mssts() const;
  void ts_ls_layer_mssts(const blitz::Array<double, 3>& ts_ls_layer_mssts_in);
  
  const blitz::Array<double, 2>& ts_ls_surf_mssts() const;
  void ts_ls_surf_mssts(const blitz::Array<double, 2>& ts_ls_surf_mssts_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linoutputs : public Lidort_Structure {
public:
  Lidort_Linoutputs();
  Lidort_Linoutputs(const Lidort_Linoutputs& src);
  ~Lidort_Linoutputs();

  const Lidort_Linatmos& atmos() const;
  void atmos(Lidort_Linatmos& atmos_in);
  
  const Lidort_Linsurf& surf() const;
  void surf(Lidort_Linsurf& surf_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Brdf : public Lidort_Structure {
public:
  Lidort_Linsup_Brdf();
  Lidort_Linsup_Brdf(const Lidort_Linsup_Brdf& src);
  ~Lidort_Linsup_Brdf();

  const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc() const;
  void ts_ls_exactdb_brdfunc(const blitz::Array<double, 4>& ts_ls_exactdb_brdfunc_in);
  
  const blitz::Array<double, 4>& ts_ls_brdf_f_0() const;
  void ts_ls_brdf_f_0(const blitz::Array<double, 4>& ts_ls_brdf_f_0_in);
  
  const blitz::Array<double, 4>& ts_ls_brdf_f() const;
  void ts_ls_brdf_f(const blitz::Array<double, 4>& ts_ls_brdf_f_in);
  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f_0() const;
  void ts_ls_user_brdf_f_0(const blitz::Array<double, 4>& ts_ls_user_brdf_f_0_in);
  
  const blitz::Array<double, 4>& ts_ls_user_brdf_f() const;
  void ts_ls_user_brdf_f(const blitz::Array<double, 4>& ts_ls_user_brdf_f_in);
  
  const blitz::Array<double, 2>& ts_ls_emissivity() const;
  void ts_ls_emissivity(const blitz::Array<double, 2>& ts_ls_emissivity_in);
  
  const blitz::Array<double, 2>& ts_ls_user_emissivity() const;
  void ts_ls_user_emissivity(const blitz::Array<double, 2>& ts_ls_user_emissivity_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Sleave : public Lidort_Structure {
public:
  Lidort_Linsup_Sleave();
  Lidort_Linsup_Sleave(const Lidort_Linsup_Sleave& src);
  ~Lidort_Linsup_Sleave();

  const blitz::Array<double, 2>& ts_lssl_slterm_isotropic() const;
  void ts_lssl_slterm_isotropic(const blitz::Array<double, 2>& ts_lssl_slterm_isotropic_in);
  
  const blitz::Array<double, 4>& ts_lssl_slterm_userangles() const;
  void ts_lssl_slterm_userangles(const blitz::Array<double, 4>& ts_lssl_slterm_userangles_in);
  
  const blitz::Array<double, 4>& ts_lssl_slterm_f_0() const;
  void ts_lssl_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_slterm_f_0_in);
  
  const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0() const;
  void ts_lssl_user_slterm_f_0(const blitz::Array<double, 4>& ts_lssl_user_slterm_f_0_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Ss_Atmos : public Lidort_Structure {
public:
  Lidort_Linsup_Ss_Atmos();
  Lidort_Linsup_Ss_Atmos(const Lidort_Linsup_Ss_Atmos& src);
  ~Lidort_Linsup_Ss_Atmos();

  const blitz::Array<double, 4>& ts_columnwf_ss() const;
  void ts_columnwf_ss(const blitz::Array<double, 4>& ts_columnwf_ss_in);
  
  const blitz::Array<double, 3>& ts_columnwf_db() const;
  void ts_columnwf_db(const blitz::Array<double, 3>& ts_columnwf_db_in);
  
  const blitz::Array<double, 5>& ts_profilewf_ss() const;
  void ts_profilewf_ss(const blitz::Array<double, 5>& ts_profilewf_ss_in);
  
  const blitz::Array<double, 4>& ts_profilewf_db() const;
  void ts_profilewf_db(const blitz::Array<double, 4>& ts_profilewf_db_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Ss_Surf : public Lidort_Structure {
public:
  Lidort_Linsup_Ss_Surf();
  Lidort_Linsup_Ss_Surf(const Lidort_Linsup_Ss_Surf& src);
  ~Lidort_Linsup_Ss_Surf();

  const blitz::Array<double, 3>& ts_surfacewf_db() const;
  void ts_surfacewf_db(const blitz::Array<double, 3>& ts_surfacewf_db_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Ss : public Lidort_Structure {
public:
  Lidort_Linsup_Ss();
  Lidort_Linsup_Ss(const Lidort_Linsup_Ss& src);
  ~Lidort_Linsup_Ss();

  const Lidort_Linsup_Ss_Atmos& atmos() const;
  void atmos(Lidort_Linsup_Ss_Atmos& atmos_in);
  
  const Lidort_Linsup_Ss_Surf& surf() const;
  void surf(Lidort_Linsup_Ss_Surf& surf_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Linsup_Inout : public Lidort_Structure {
public:
  Lidort_Linsup_Inout();
  Lidort_Linsup_Inout(const Lidort_Linsup_Inout& src);
  ~Lidort_Linsup_Inout();

  const Lidort_Linsup_Brdf& brdf() const;
  void brdf(Lidort_Linsup_Brdf& brdf_in);
  
  const Lidort_Linsup_Ss& ss() const;
  void ss(Lidort_Linsup_Ss& ss_in);
  
  const Lidort_Linsup_Sleave& sleave() const;
  void sleave(Lidort_Linsup_Sleave& sleave_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Main_Outputs : public Lidort_Structure {
public:
  Lidort_Main_Outputs();
  Lidort_Main_Outputs(const Lidort_Main_Outputs& src);
  ~Lidort_Main_Outputs();

  const blitz::Array<double, 3>& ts_intensity() const;
  void ts_intensity(const blitz::Array<double, 3>& ts_intensity_in);
  
  const blitz::Array<double, 3>& ts_meani_diffuse() const;
  void ts_meani_diffuse(const blitz::Array<double, 3>& ts_meani_diffuse_in);
  
  const blitz::Array<double, 3>& ts_flux_diffuse() const;
  void ts_flux_diffuse(const blitz::Array<double, 3>& ts_flux_diffuse_in);
  
  const blitz::Array<double, 2>& ts_dnmeani_direct() const;
  void ts_dnmeani_direct(const blitz::Array<double, 2>& ts_dnmeani_direct_in);
  
  const blitz::Array<double, 2>& ts_dnflux_direct() const;
  void ts_dnflux_direct(const blitz::Array<double, 2>& ts_dnflux_direct_in);
  
  const blitz::Array<double, 1>& ts_albmed_user() const;
  void ts_albmed_user(const blitz::Array<double, 1>& ts_albmed_user_in);
  
  const blitz::Array<double, 1>& ts_trnmed_user() const;
  void ts_trnmed_user(const blitz::Array<double, 1>& ts_trnmed_user_in);
  
  const blitz::Array<double, 1>& ts_albmed_fluxes() const;
  void ts_albmed_fluxes(const blitz::Array<double, 1>& ts_albmed_fluxes_in);
  
  const blitz::Array<double, 1>& ts_trnmed_fluxes() const;
  void ts_trnmed_fluxes(const blitz::Array<double, 1>& ts_trnmed_fluxes_in);
  
  const blitz::Array<double, 1>& ts_planetary_transterm() const;
  void ts_planetary_transterm(const blitz::Array<double, 1>& ts_planetary_transterm_in);
  
  const double& ts_planetary_sbterm() const;
  void ts_planetary_sbterm(const double& ts_planetary_sbterm_in);
  
  const blitz::Array<double, 2>& ts_pathgeoms() const;
  void ts_pathgeoms(const blitz::Array<double, 2>& ts_pathgeoms_in);
  
  const blitz::Array<double, 2>& ts_lostrans() const;
  void ts_lostrans(const blitz::Array<double, 2>& ts_lostrans_in);
  
  const blitz::Array<double, 2>& ts_layer_mssts() const;
  void ts_layer_mssts(const blitz::Array<double, 2>& ts_layer_mssts_in);
  
  const blitz::Array<double, 1>& ts_surf_mssts() const;
  void ts_surf_mssts(const blitz::Array<double, 1>& ts_surf_mssts_in);
  
  const blitz::Array<double, 2>& ts_contribs() const;
  void ts_contribs(const blitz::Array<double, 2>& ts_contribs_in);
  
  const blitz::Array<int, 1>& ts_fourier_saved() const;
  void ts_fourier_saved(const blitz::Array<int, 1>& ts_fourier_saved_in);
  
  const int& ts_n_geometries() const;
  void ts_n_geometries(const int& ts_n_geometries_in);
  
  const blitz::Array<double, 1>& ts_solarbeam_boatrans() const;
  void ts_solarbeam_boatrans(const blitz::Array<double, 1>& ts_solarbeam_boatrans_in);
  
  const double& ts_spheralb() const;
  void ts_spheralb(const double& ts_spheralb_in);
  
  const blitz::Array<double, 1>& ts_trans1_user() const;
  void ts_trans1_user(const blitz::Array<double, 1>& ts_trans1_user_in);
  
  const blitz::Array<double, 1>& ts_trans1_beam() const;
  void ts_trans1_beam(const blitz::Array<double, 1>& ts_trans1_beam_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Wladjusted_Outputs : public Lidort_Structure {
public:
  Lidort_Wladjusted_Outputs();
  Lidort_Wladjusted_Outputs(const Lidort_Wladjusted_Outputs& src);
  ~Lidort_Wladjusted_Outputs();

  const blitz::Array<double, 1>& ts_wladjusted_isotropic() const;
  void ts_wladjusted_isotropic(const blitz::Array<double, 1>& ts_wladjusted_isotropic_in);
  
  const blitz::Array<double, 3>& ts_wladjusted_direct() const;
  void ts_wladjusted_direct(const blitz::Array<double, 3>& ts_wladjusted_direct_in);
  
  const blitz::Array<double, 3>& ts_wladjusted_f_ords_0() const;
  void ts_wladjusted_f_ords_0(const blitz::Array<double, 3>& ts_wladjusted_f_ords_0_in);
  
  const blitz::Array<double, 3>& ts_wladjusted_f_user_0() const;
  void ts_wladjusted_f_user_0(const blitz::Array<double, 3>& ts_wladjusted_f_user_0_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Exception_Handling : public Lidort_Structure {
public:
  Lidort_Exception_Handling();
  Lidort_Exception_Handling(const Lidort_Exception_Handling& src);
  ~Lidort_Exception_Handling();

  const int& ts_status_inputcheck() const;
  void ts_status_inputcheck(const int& ts_status_inputcheck_in);
  
  const int& ts_ncheckmessages() const;
  void ts_ncheckmessages(const int& ts_ncheckmessages_in);
  
  const std::vector< std::string > ts_checkmessages() const;
  
  const std::vector< std::string > ts_actions() const;
  
  const int& ts_status_calculation() const;
  void ts_status_calculation(const int& ts_status_calculation_in);
  
  const std::string ts_message() const;
  
  const std::string ts_trace_1() const;
  
  const std::string ts_trace_2() const;
  
  const std::string ts_trace_3() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Input_Exception_Handling : public Lidort_Structure {
public:
  Lidort_Input_Exception_Handling();
  Lidort_Input_Exception_Handling(const Lidort_Input_Exception_Handling& src);
  ~Lidort_Input_Exception_Handling();

  const int& ts_status_inputread() const;
  void ts_status_inputread(const int& ts_status_inputread_in);
  
  const int& ts_ninputmessages() const;
  void ts_ninputmessages(const int& ts_ninputmessages_in);
  
  const std::vector< std::string > ts_inputmessages() const;
  
  const std::vector< std::string > ts_inputactions() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Outputs : public Lidort_Structure {
public:
  Lidort_Outputs();
  Lidort_Outputs(const Lidort_Outputs& src);
  ~Lidort_Outputs();

  const Lidort_Main_Outputs& main() const;
  void main(Lidort_Main_Outputs& main_in);
  
  const Lidort_Wladjusted_Outputs& wlout() const;
  void wlout(Lidort_Wladjusted_Outputs& wlout_in);
  
  const Lidort_Exception_Handling& status() const;
  void status(Lidort_Exception_Handling& status_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Sup_Brdf : public Lidort_Structure {
public:
  Lidort_Sup_Brdf();
  Lidort_Sup_Brdf(const Lidort_Sup_Brdf& src);
  ~Lidort_Sup_Brdf();

  const blitz::Array<double, 3>& ts_exactdb_brdfunc() const;
  void ts_exactdb_brdfunc(const blitz::Array<double, 3>& ts_exactdb_brdfunc_in);
  
  const blitz::Array<double, 3>& ts_brdf_f_0() const;
  void ts_brdf_f_0(const blitz::Array<double, 3>& ts_brdf_f_0_in);
  
  const blitz::Array<double, 3>& ts_brdf_f() const;
  void ts_brdf_f(const blitz::Array<double, 3>& ts_brdf_f_in);
  
  const blitz::Array<double, 3>& ts_user_brdf_f_0() const;
  void ts_user_brdf_f_0(const blitz::Array<double, 3>& ts_user_brdf_f_0_in);
  
  const blitz::Array<double, 3>& ts_user_brdf_f() const;
  void ts_user_brdf_f(const blitz::Array<double, 3>& ts_user_brdf_f_in);
  
  const blitz::Array<double, 1>& ts_emissivity() const;
  void ts_emissivity(const blitz::Array<double, 1>& ts_emissivity_in);
  
  const blitz::Array<double, 1>& ts_user_emissivity() const;
  void ts_user_emissivity(const blitz::Array<double, 1>& ts_user_emissivity_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Sup_Sleave : public Lidort_Structure {
public:
  Lidort_Sup_Sleave();
  Lidort_Sup_Sleave(const Lidort_Sup_Sleave& src);
  ~Lidort_Sup_Sleave();

  const blitz::Array<double, 1>& ts_slterm_isotropic() const;
  void ts_slterm_isotropic(const blitz::Array<double, 1>& ts_slterm_isotropic_in);
  
  const blitz::Array<double, 3>& ts_slterm_userangles() const;
  void ts_slterm_userangles(const blitz::Array<double, 3>& ts_slterm_userangles_in);
  
  const blitz::Array<double, 3>& ts_slterm_f_0() const;
  void ts_slterm_f_0(const blitz::Array<double, 3>& ts_slterm_f_0_in);
  
  const blitz::Array<double, 3>& ts_user_slterm_f_0() const;
  void ts_user_slterm_f_0(const blitz::Array<double, 3>& ts_user_slterm_f_0_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Sup_Ss : public Lidort_Structure {
public:
  Lidort_Sup_Ss();
  Lidort_Sup_Ss(const Lidort_Sup_Ss& src);
  ~Lidort_Sup_Ss();

  const blitz::Array<double, 3>& ts_intensity_ss() const;
  void ts_intensity_ss(const blitz::Array<double, 3>& ts_intensity_ss_in);
  
  const blitz::Array<double, 2>& ts_intensity_db() const;
  void ts_intensity_db(const blitz::Array<double, 2>& ts_intensity_db_in);
  
  const blitz::Array<double, 2>& ts_contribs_ss() const;
  void ts_contribs_ss(const blitz::Array<double, 2>& ts_contribs_ss_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Sup_Inout : public Lidort_Structure {
public:
  Lidort_Sup_Inout();
  Lidort_Sup_Inout(const Lidort_Sup_Inout& src);
  ~Lidort_Sup_Inout();

  const Lidort_Sup_Brdf& brdf() const;
  void brdf(Lidort_Sup_Brdf& brdf_in);
  
  const Lidort_Sup_Ss& ss() const;
  void ss(Lidort_Sup_Ss& ss_in);
  
  const Lidort_Sup_Sleave& sleave() const;
  void sleave(Lidort_Sup_Sleave& sleave_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Boolean : public Lidort_Structure {
public:
  Lidort_Fixed_Boolean();
  Lidort_Fixed_Boolean(const Lidort_Fixed_Boolean& src);
  ~Lidort_Fixed_Boolean();

  const bool ts_do_fullrad_mode() const;
  void ts_do_fullrad_mode(const bool& ts_do_fullrad_mode_in);
  
  const bool ts_do_thermal_emission() const;
  void ts_do_thermal_emission(const bool& ts_do_thermal_emission_in);
  
  const bool ts_do_surface_emission() const;
  void ts_do_surface_emission(const bool& ts_do_surface_emission_in);
  
  const bool ts_do_plane_parallel() const;
  void ts_do_plane_parallel(const bool& ts_do_plane_parallel_in);
  
  const bool ts_do_brdf_surface() const;
  void ts_do_brdf_surface(const bool& ts_do_brdf_surface_in);
  
  const bool ts_do_upwelling() const;
  void ts_do_upwelling(const bool& ts_do_upwelling_in);
  
  const bool ts_do_dnwelling() const;
  void ts_do_dnwelling(const bool& ts_do_dnwelling_in);
  
  const bool ts_do_toa_contribs() const;
  void ts_do_toa_contribs(const bool& ts_do_toa_contribs_in);
  
  const bool ts_do_surface_leaving() const;
  void ts_do_surface_leaving(const bool& ts_do_surface_leaving_in);
  
  const bool ts_do_sl_isotropic() const;
  void ts_do_sl_isotropic(const bool& ts_do_sl_isotropic_in);
  
  const bool ts_do_water_leaving() const;
  void ts_do_water_leaving(const bool& ts_do_water_leaving_in);
  
  const bool ts_do_fluorescence() const;
  void ts_do_fluorescence(const bool& ts_do_fluorescence_in);
  
  const bool ts_do_tf_iteration() const;
  void ts_do_tf_iteration(const bool& ts_do_tf_iteration_in);
  
  const bool ts_do_wladjusted_output() const;
  void ts_do_wladjusted_output(const bool& ts_do_wladjusted_output_in);
  
  const bool ts_do_toa_illumination() const;
  void ts_do_toa_illumination(const bool& ts_do_toa_illumination_in);
  
  const bool ts_do_boa_illumination() const;
  void ts_do_boa_illumination(const bool& ts_do_boa_illumination_in);
  
  const blitz::Array<bool, 1> ts_do_albtrn_media() const;
  void ts_do_albtrn_media(const blitz::Array<bool, 1>& ts_do_albtrn_media_in);
  
  const bool ts_do_planetary_problem() const;
  void ts_do_planetary_problem(const bool& ts_do_planetary_problem_in);
  
  const bool ts_do_mssts() const;
  void ts_do_mssts(const bool& ts_do_mssts_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Control : public Lidort_Structure {
public:
  Lidort_Fixed_Control();
  Lidort_Fixed_Control(const Lidort_Fixed_Control& src);
  ~Lidort_Fixed_Control();

  const int& ts_taylor_order() const;
  void ts_taylor_order(const int& ts_taylor_order_in);
  
  const int& ts_nstreams() const;
  void ts_nstreams(const int& ts_nstreams_in);
  
  const int& ts_nlayers() const;
  void ts_nlayers(const int& ts_nlayers_in);
  
  const int& ts_nfinelayers() const;
  void ts_nfinelayers(const int& ts_nfinelayers_in);
  
  const int& ts_n_thermal_coeffs() const;
  void ts_n_thermal_coeffs(const int& ts_n_thermal_coeffs_in);
  
  const double& ts_lidort_accuracy() const;
  void ts_lidort_accuracy(const double& ts_lidort_accuracy_in);
  
  const double& ts_asymtx_tolerance() const;
  void ts_asymtx_tolerance(const double& ts_asymtx_tolerance_in);
  
  const int& ts_tf_maxiter() const;
  void ts_tf_maxiter(const int& ts_tf_maxiter_in);
  
  const double& ts_tf_criterion() const;
  void ts_tf_criterion(const double& ts_tf_criterion_in);
  
  const double& ts_toa_illumination() const;
  void ts_toa_illumination(const double& ts_toa_illumination_in);
  
  const double& ts_boa_illumination() const;
  void ts_boa_illumination(const double& ts_boa_illumination_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Sunrays : public Lidort_Structure {
public:
  Lidort_Fixed_Sunrays();
  Lidort_Fixed_Sunrays(const Lidort_Fixed_Sunrays& src);
  ~Lidort_Fixed_Sunrays();

  const double& ts_flux_factor() const;
  void ts_flux_factor(const double& ts_flux_factor_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Uservalues : public Lidort_Structure {
public:
  Lidort_Fixed_Uservalues();
  Lidort_Fixed_Uservalues(const Lidort_Fixed_Uservalues& src);
  ~Lidort_Fixed_Uservalues();

  const int& ts_n_user_levels() const;
  void ts_n_user_levels(const int& ts_n_user_levels_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Chapman : public Lidort_Structure {
public:
  Lidort_Fixed_Chapman();
  Lidort_Fixed_Chapman(const Lidort_Fixed_Chapman& src);
  ~Lidort_Fixed_Chapman();

  const blitz::Array<double, 1>& ts_height_grid() const;
  void ts_height_grid(const blitz::Array<double, 1>& ts_height_grid_in);
  
  const blitz::Array<double, 1>& ts_pressure_grid() const;
  void ts_pressure_grid(const blitz::Array<double, 1>& ts_pressure_grid_in);
  
  const blitz::Array<double, 1>& ts_temperature_grid() const;
  void ts_temperature_grid(const blitz::Array<double, 1>& ts_temperature_grid_in);
  
  const blitz::Array<int, 1>& ts_finegrid() const;
  void ts_finegrid(const blitz::Array<int, 1>& ts_finegrid_in);
  
  const double& ts_rfindex_parameter() const;
  void ts_rfindex_parameter(const double& ts_rfindex_parameter_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Optical : public Lidort_Structure {
public:
  Lidort_Fixed_Optical();
  Lidort_Fixed_Optical(const Lidort_Fixed_Optical& src);
  ~Lidort_Fixed_Optical();

  const blitz::Array<double, 1>& ts_deltau_vert_input() const;
  void ts_deltau_vert_input(const blitz::Array<double, 1>& ts_deltau_vert_input_in);
  
  const blitz::Array<double, 2>& ts_phasmoms_total_input() const;
  void ts_phasmoms_total_input(const blitz::Array<double, 2>& ts_phasmoms_total_input_in);
  
  const blitz::Array<double, 2>& ts_phasfunc_input_up() const;
  void ts_phasfunc_input_up(const blitz::Array<double, 2>& ts_phasfunc_input_up_in);
  
  const blitz::Array<double, 2>& ts_phasfunc_input_dn() const;
  void ts_phasfunc_input_dn(const blitz::Array<double, 2>& ts_phasfunc_input_dn_in);
  
  const double& ts_lambertian_albedo() const;
  void ts_lambertian_albedo(const double& ts_lambertian_albedo_in);
  
  const blitz::Array<double, 1>& ts_thermal_bb_input() const;
  void ts_thermal_bb_input(const blitz::Array<double, 1>& ts_thermal_bb_input_in);
  
  const double& ts_surface_bb_input() const;
  void ts_surface_bb_input(const double& ts_surface_bb_input_in);
  
  const double& ts_atmos_wavelength() const;
  void ts_atmos_wavelength(const double& ts_atmos_wavelength_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Write : public Lidort_Structure {
public:
  Lidort_Fixed_Write();
  Lidort_Fixed_Write(const Lidort_Fixed_Write& src);
  ~Lidort_Fixed_Write();

  const bool ts_do_debug_write() const;
  void ts_do_debug_write(const bool& ts_do_debug_write_in);
  
  const bool ts_do_write_input() const;
  void ts_do_write_input(const bool& ts_do_write_input_in);
  
  const std::string ts_input_write_filename() const;
  
  const bool ts_do_write_scenario() const;
  void ts_do_write_scenario(const bool& ts_do_write_scenario_in);
  
  const std::string ts_scenario_write_filename() const;
  
  const bool ts_do_write_fourier() const;
  void ts_do_write_fourier(const bool& ts_do_write_fourier_in);
  
  const std::string ts_fourier_write_filename() const;
  
  const bool ts_do_write_results() const;
  void ts_do_write_results(const bool& ts_do_write_results_in);
  
  const std::string ts_results_write_filename() const;
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Fixed_Inputs : public Lidort_Structure {
public:
  Lidort_Fixed_Inputs();
  Lidort_Fixed_Inputs(const Lidort_Fixed_Inputs& src);
  ~Lidort_Fixed_Inputs();

  const Lidort_Fixed_Boolean& f_bool() const;
  void f_bool(Lidort_Fixed_Boolean& bool_in);
  
  const Lidort_Fixed_Control& cont() const;
  void cont(Lidort_Fixed_Control& cont_in);
  
  const Lidort_Fixed_Sunrays& sunrays() const;
  void sunrays(Lidort_Fixed_Sunrays& sunrays_in);
  
  const Lidort_Fixed_Uservalues& userval() const;
  void userval(Lidort_Fixed_Uservalues& userval_in);
  
  const Lidort_Fixed_Chapman& chapman() const;
  void chapman(Lidort_Fixed_Chapman& chapman_in);
  
  const Lidort_Fixed_Optical& optical() const;
  void optical(Lidort_Fixed_Optical& optical_in);
  
  const Lidort_Fixed_Write& write() const;
  void write(Lidort_Fixed_Write& write_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Boolean : public Lidort_Structure {
public:
  Lidort_Modified_Boolean();
  Lidort_Modified_Boolean(const Lidort_Modified_Boolean& src);
  ~Lidort_Modified_Boolean();

  const bool ts_do_focorr() const;
  void ts_do_focorr(const bool& ts_do_focorr_in);
  
  const bool ts_do_focorr_external() const;
  void ts_do_focorr_external(const bool& ts_do_focorr_external_in);
  
  const bool ts_do_focorr_nadir() const;
  void ts_do_focorr_nadir(const bool& ts_do_focorr_nadir_in);
  
  const bool ts_do_focorr_outgoing() const;
  void ts_do_focorr_outgoing(const bool& ts_do_focorr_outgoing_in);
  
  const bool ts_do_sscorr_truncation() const;
  void ts_do_sscorr_truncation(const bool& ts_do_sscorr_truncation_in);
  
  const bool ts_do_sscorr_usephasfunc() const;
  void ts_do_sscorr_usephasfunc(const bool& ts_do_sscorr_usephasfunc_in);
  
  const bool ts_do_external_wleave() const;
  void ts_do_external_wleave(const bool& ts_do_external_wleave_in);
  
  const bool ts_do_double_convtest() const;
  void ts_do_double_convtest(const bool& ts_do_double_convtest_in);
  
  const bool ts_do_solar_sources() const;
  void ts_do_solar_sources(const bool& ts_do_solar_sources_in);
  
  const bool ts_do_refractive_geometry() const;
  void ts_do_refractive_geometry(const bool& ts_do_refractive_geometry_in);
  
  const bool ts_do_chapman_function() const;
  void ts_do_chapman_function(const bool& ts_do_chapman_function_in);
  
  const bool ts_do_rayleigh_only() const;
  void ts_do_rayleigh_only(const bool& ts_do_rayleigh_only_in);
  
  const bool ts_do_isotropic_only() const;
  void ts_do_isotropic_only(const bool& ts_do_isotropic_only_in);
  
  const bool ts_do_no_azimuth() const;
  void ts_do_no_azimuth(const bool& ts_do_no_azimuth_in);
  
  const bool ts_do_all_fourier() const;
  void ts_do_all_fourier(const bool& ts_do_all_fourier_in);
  
  const bool ts_do_deltam_scaling() const;
  void ts_do_deltam_scaling(const bool& ts_do_deltam_scaling_in);
  
  const bool ts_do_solution_saving() const;
  void ts_do_solution_saving(const bool& ts_do_solution_saving_in);
  
  const bool ts_do_bvp_telescoping() const;
  void ts_do_bvp_telescoping(const bool& ts_do_bvp_telescoping_in);
  
  const bool ts_do_user_streams() const;
  void ts_do_user_streams(const bool& ts_do_user_streams_in);
  
  const bool ts_do_additional_mvout() const;
  void ts_do_additional_mvout(const bool& ts_do_additional_mvout_in);
  
  const bool ts_do_mvout_only() const;
  void ts_do_mvout_only(const bool& ts_do_mvout_only_in);
  
  const bool ts_do_thermal_transonly() const;
  void ts_do_thermal_transonly(const bool& ts_do_thermal_transonly_in);
  
  const bool ts_do_observation_geometry() const;
  void ts_do_observation_geometry(const bool& ts_do_observation_geometry_in);
  
  const bool ts_do_doublet_geometry() const;
  void ts_do_doublet_geometry(const bool& ts_do_doublet_geometry_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Control : public Lidort_Structure {
public:
  Lidort_Modified_Control();
  Lidort_Modified_Control(const Lidort_Modified_Control& src);
  ~Lidort_Modified_Control();

  const int& ts_nmoments_input() const;
  void ts_nmoments_input(const int& ts_nmoments_input_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Sunrays : public Lidort_Structure {
public:
  Lidort_Modified_Sunrays();
  Lidort_Modified_Sunrays(const Lidort_Modified_Sunrays& src);
  ~Lidort_Modified_Sunrays();

  const int& ts_nbeams() const;
  void ts_nbeams(const int& ts_nbeams_in);
  
  const blitz::Array<double, 1>& ts_beam_szas() const;
  void ts_beam_szas(const blitz::Array<double, 1>& ts_beam_szas_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Uservalues : public Lidort_Structure {
public:
  Lidort_Modified_Uservalues();
  Lidort_Modified_Uservalues(const Lidort_Modified_Uservalues& src);
  ~Lidort_Modified_Uservalues();

  const int& ts_n_user_relazms() const;
  void ts_n_user_relazms(const int& ts_n_user_relazms_in);
  
  const blitz::Array<double, 1>& ts_user_relazms() const;
  void ts_user_relazms(const blitz::Array<double, 1>& ts_user_relazms_in);
  
  const int& ts_n_user_streams() const;
  void ts_n_user_streams(const int& ts_n_user_streams_in);
  
  const blitz::Array<double, 1>& ts_user_angles_input() const;
  void ts_user_angles_input(const blitz::Array<double, 1>& ts_user_angles_input_in);
  
  const blitz::Array<double, 1>& ts_user_levels() const;
  void ts_user_levels(const blitz::Array<double, 1>& ts_user_levels_in);
  
  const double& ts_geometry_specheight() const;
  void ts_geometry_specheight(const double& ts_geometry_specheight_in);
  
  const int& ts_n_user_obsgeoms() const;
  void ts_n_user_obsgeoms(const int& ts_n_user_obsgeoms_in);
  
  const blitz::Array<double, 2>& ts_user_obsgeoms_input() const;
  void ts_user_obsgeoms_input(const blitz::Array<double, 2>& ts_user_obsgeoms_input_in);
  
  const int& ts_n_user_doublets() const;
  void ts_n_user_doublets(const int& ts_n_user_doublets_in);
  
  const blitz::Array<double, 2>& ts_user_doublets() const;
  void ts_user_doublets(const blitz::Array<double, 2>& ts_user_doublets_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Chapman : public Lidort_Structure {
public:
  Lidort_Modified_Chapman();
  Lidort_Modified_Chapman(const Lidort_Modified_Chapman& src);
  ~Lidort_Modified_Chapman();

  const double& ts_earth_radius() const;
  void ts_earth_radius(const double& ts_earth_radius_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Optical : public Lidort_Structure {
public:
  Lidort_Modified_Optical();
  Lidort_Modified_Optical(const Lidort_Modified_Optical& src);
  ~Lidort_Modified_Optical();

  const blitz::Array<double, 1>& ts_omega_total_input() const;
  void ts_omega_total_input(const blitz::Array<double, 1>& ts_omega_total_input_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};

class Lidort_Modified_Inputs : public Lidort_Structure {
public:
  Lidort_Modified_Inputs();
  Lidort_Modified_Inputs(const Lidort_Modified_Inputs& src);
  ~Lidort_Modified_Inputs();

  const Lidort_Modified_Boolean& mbool() const;
  void mbool(Lidort_Modified_Boolean& mbool_in);
  
  const Lidort_Modified_Control& mcont() const;
  void mcont(Lidort_Modified_Control& mcont_in);
  
  const Lidort_Modified_Sunrays& msunrays() const;
  void msunrays(Lidort_Modified_Sunrays& msunrays_in);
  
  const Lidort_Modified_Uservalues& muserval() const;
  void muserval(Lidort_Modified_Uservalues& muserval_in);
  
  const Lidort_Modified_Chapman& mchapman() const;
  void mchapman(Lidort_Modified_Chapman& mchapman_in);
  
  const Lidort_Modified_Optical& moptical() const;
  void moptical(Lidort_Modified_Optical& moptical_in);
  
  

  virtual void print(std::ostream &output_stream) const;
};



}
