module lidort_interface_types

use iso_c_binding
use lidort_pars_m

! This module was auto-generated 

implicit none


! Links to module: "lidort_pars_m" in file: "lidort_pars.F90.in"
type, bind(c) :: lidort_pars_m_c
  character(kind=c_char) :: lidort_version_number(5)
  integer(c_int) :: lidort_inunit
  integer(c_int) :: lidort_scenunit
  integer(c_int) :: lidort_funit
  integer(c_int) :: lidort_resunit
  integer(c_int) :: lidort_errunit
  integer(c_int) :: lidort_dbgunit
  integer(c_int) :: max_messages
  integer(c_int) :: maxstreams
  integer(c_int) :: maxlayers
  integer(c_int) :: maxfinelayers
  integer(c_int) :: maxmoments_input
  integer(c_int) :: max_thermal_coeffs
  integer(c_int) :: maxbeams
  integer(c_int) :: max_user_streams
  integer(c_int) :: max_user_relazms
  integer(c_int) :: max_user_obsgeoms
  integer(c_int) :: max_user_levels
  integer(c_int) :: max_partlayers
  integer(c_int) :: max_taylor_terms
  integer(c_int) :: max_directions
  integer(c_int) :: max_brdf_kernels
  integer(c_int) :: max_brdf_parameters
  integer(c_int) :: maxstreams_brdf
  integer(c_int) :: max_msrs_muquad
  integer(c_int) :: max_msrs_phiquad
  integer(c_int) :: maxstreams_scaling
  integer(c_int) :: max_atmoswfs
  integer(c_int) :: max_surfacewfs
  integer(c_int) :: max_sleavewfs
  integer(c_int) :: max_geometries
  integer(c_int) :: max_allstrms
  integer(c_int) :: max_allstrms_p1
  integer(c_int) :: maxmoments
  integer(c_int) :: maxfourier
  integer(c_int) :: maxsthalf_brdf
  integer(c_int) :: maxstreams_2
  integer(c_int) :: maxstreams_p1
  integer(c_int) :: maxtotal
  integer(c_int) :: maxbandtotal
  real(c_double) :: one
  real(c_double) :: zero
  real(c_double) :: onep5
  real(c_double) :: two
  real(c_double) :: three
  real(c_double) :: four
  real(c_double) :: quarter
  real(c_double) :: half
  real(c_double) :: minus_one
  real(c_double) :: minus_two
  real(c_double) :: pie
  real(c_double) :: deg_to_rad
  real(c_double) :: pi2
  real(c_double) :: pi4
  real(c_double) :: pio2
  real(c_double) :: pio4
  real(c_double) :: eps3
  real(c_double) :: eps4
  real(c_double) :: eps5
  real(c_double) :: taylor_small
  real(c_double) :: smallnum
  real(c_double) :: bigexp
  real(c_double) :: hopital_tolerance
  real(c_double) :: omega_smallnum
  real(c_double) :: max_tau_spath
  real(c_double) :: max_tau_upath
  real(c_double) :: max_tau_qpath
  integer(c_int) :: lidort_serious
  integer(c_int) :: lidort_warning
  integer(c_int) :: lidort_info
  integer(c_int) :: lidort_debug
  integer(c_int) :: lidort_success
  integer(c_int) :: upidx
  integer(c_int) :: dnidx
  integer(c_int) :: lambertian_idx
  integer(c_int) :: rossthin_idx
  integer(c_int) :: rossthick_idx
  integer(c_int) :: lisparse_idx
  integer(c_int) :: lidense_idx
  integer(c_int) :: hapke_idx
  integer(c_int) :: roujean_idx
  integer(c_int) :: rahman_idx
  integer(c_int) :: coxmunk_idx
  integer(c_int) :: bpdfsoil_idx
  integer(c_int) :: bpdfvegn_idx
  integer(c_int) :: bpdfndvi_idx
  integer(c_int) :: newcmglint_idx
  integer(c_int) :: rtkhotspot_idx
  integer(c_int) :: modfresnel_idx
  integer(c_int) :: snowbrdf_idx
  integer(c_int) :: maxbrdf_idx
  
end type lidort_pars_m_c



! Links to type: "brdf_linsup_inputs" from module: "brdf_lin_sup_inputs_def_m" in file: "brdf_lin_sup_inputs_def.F90"
type, bind(c) :: brdf_linsup_inputs_c
  type(c_ptr) :: bs_do_kernel_factor_wfs ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_do_kernel_factor_wfs_f_shapes
  integer(c_int) :: bs_do_kernel_factor_wfs_f_byte_size

  type(c_ptr) :: bs_do_kernel_params_wfs ! dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS)
  integer(c_int), dimension(2) :: bs_do_kernel_params_wfs_f_shapes
  integer(c_int) :: bs_do_kernel_params_wfs_f_byte_size

  type(c_ptr) :: bs_do_kparams_derivs ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_do_kparams_derivs_f_shapes
  integer(c_int) :: bs_do_kparams_derivs_f_byte_size

  type(c_ptr) :: bs_n_surface_wfs ! scalar
  integer(c_int) :: bs_n_surface_wfs_f_byte_size

  type(c_ptr) :: bs_n_kernel_factor_wfs ! scalar
  integer(c_int) :: bs_n_kernel_factor_wfs_f_byte_size

  type(c_ptr) :: bs_n_kernel_params_wfs ! scalar
  integer(c_int) :: bs_n_kernel_params_wfs_f_byte_size

  type(c_ptr) :: bs_do_bsavalue_wf ! scalar
  integer(c_int) :: bs_do_bsavalue_wf_f_byte_size

  type(c_ptr) :: bs_do_wsavalue_wf ! scalar
  integer(c_int) :: bs_do_wsavalue_wf_f_byte_size

  type(c_ptr) :: bs_do_windspeed_wf ! scalar
  integer(c_int) :: bs_do_windspeed_wf_f_byte_size

  
end type brdf_linsup_inputs_c

! Links to type: "brdf_linsup_outputs" from module: "brdf_lin_sup_outputs_def_m" in file: "brdf_lin_sup_outputs_def.F90"
type, bind(c) :: brdf_linsup_outputs_c
  type(c_ptr) :: bs_ls_dbounce_brdfunc ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(4) :: bs_ls_dbounce_brdfunc_f_shapes
  integer(c_int) :: bs_ls_dbounce_brdfunc_f_byte_size

  type(c_ptr) :: bs_ls_brdf_f_0 ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: bs_ls_brdf_f_0_f_shapes
  integer(c_int) :: bs_ls_brdf_f_0_f_byte_size

  type(c_ptr) :: bs_ls_brdf_f ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS)
  integer(c_int), dimension(4) :: bs_ls_brdf_f_f_shapes
  integer(c_int) :: bs_ls_brdf_f_f_byte_size

  type(c_ptr) :: bs_ls_user_brdf_f_0 ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: bs_ls_user_brdf_f_0_f_shapes
  integer(c_int) :: bs_ls_user_brdf_f_0_f_byte_size

  type(c_ptr) :: bs_ls_user_brdf_f ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS)
  integer(c_int), dimension(4) :: bs_ls_user_brdf_f_f_shapes
  integer(c_int) :: bs_ls_user_brdf_f_f_byte_size

  type(c_ptr) :: bs_ls_emissivity ! dimension(MAX_SURFACEWFS, MAXSTREAMS)
  integer(c_int), dimension(2) :: bs_ls_emissivity_f_shapes
  integer(c_int) :: bs_ls_emissivity_f_byte_size

  type(c_ptr) :: bs_ls_user_emissivity ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS)
  integer(c_int), dimension(2) :: bs_ls_user_emissivity_f_shapes
  integer(c_int) :: bs_ls_user_emissivity_f_byte_size

  
end type brdf_linsup_outputs_c

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.F90"
type, bind(c) :: brdf_sup_inputs_c
  type(c_ptr) :: bs_do_brdf_surface ! scalar
  integer(c_int) :: bs_do_brdf_surface_f_byte_size

  type(c_ptr) :: bs_do_surface_emission ! scalar
  integer(c_int) :: bs_do_surface_emission_f_byte_size

  type(c_ptr) :: bs_do_solar_sources ! scalar
  integer(c_int) :: bs_do_solar_sources_f_byte_size

  type(c_ptr) :: bs_do_user_streams ! scalar
  integer(c_int) :: bs_do_user_streams_f_byte_size

  type(c_ptr) :: bs_do_user_obsgeoms ! scalar
  integer(c_int) :: bs_do_user_obsgeoms_f_byte_size

  type(c_ptr) :: bs_do_doublet_geometry ! scalar
  integer(c_int) :: bs_do_doublet_geometry_f_byte_size

  type(c_ptr) :: bs_nstreams ! scalar
  integer(c_int) :: bs_nstreams_f_byte_size

  type(c_ptr) :: bs_nbeams ! scalar
  integer(c_int) :: bs_nbeams_f_byte_size

  type(c_ptr) :: bs_beam_szas ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: bs_beam_szas_f_shapes
  integer(c_int) :: bs_beam_szas_f_byte_size

  type(c_ptr) :: bs_n_user_relazms ! scalar
  integer(c_int) :: bs_n_user_relazms_f_byte_size

  type(c_ptr) :: bs_user_relazms ! dimension(MAX_USER_RELAZMS)
  integer(c_int), dimension(1) :: bs_user_relazms_f_shapes
  integer(c_int) :: bs_user_relazms_f_byte_size

  type(c_ptr) :: bs_n_user_streams ! scalar
  integer(c_int) :: bs_n_user_streams_f_byte_size

  type(c_ptr) :: bs_user_angles_input ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: bs_user_angles_input_f_shapes
  integer(c_int) :: bs_user_angles_input_f_byte_size

  type(c_ptr) :: bs_n_user_obsgeoms ! scalar
  integer(c_int) :: bs_n_user_obsgeoms_f_byte_size

  type(c_ptr) :: bs_user_obsgeoms ! dimension(MAX_USER_OBSGEOMS, 3)
  integer(c_int), dimension(2) :: bs_user_obsgeoms_f_shapes
  integer(c_int) :: bs_user_obsgeoms_f_byte_size

  type(c_ptr) :: bs_n_user_doublets ! scalar
  integer(c_int) :: bs_n_user_doublets_f_byte_size

  type(c_ptr) :: bs_user_doublets ! dimension(MAX_USER_STREAMS, 2)
  integer(c_int), dimension(2) :: bs_user_doublets_f_shapes
  integer(c_int) :: bs_user_doublets_f_byte_size

  type(c_ptr) :: bs_n_brdf_kernels ! scalar
  integer(c_int) :: bs_n_brdf_kernels_f_byte_size

  
  integer(c_int), dimension(1) :: bs_brdf_names_f_shapes
  integer(c_int) :: bs_brdf_names_f_len

  type(c_ptr) :: bs_which_brdf ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_which_brdf_f_shapes
  integer(c_int) :: bs_which_brdf_f_byte_size

  type(c_ptr) :: bs_n_brdf_parameters ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_n_brdf_parameters_f_shapes
  integer(c_int) :: bs_n_brdf_parameters_f_byte_size

  type(c_ptr) :: bs_brdf_parameters ! dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS)
  integer(c_int), dimension(2) :: bs_brdf_parameters_f_shapes
  integer(c_int) :: bs_brdf_parameters_f_byte_size

  type(c_ptr) :: bs_lambertian_kernel_flag ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_lambertian_kernel_flag_f_shapes
  integer(c_int) :: bs_lambertian_kernel_flag_f_byte_size

  type(c_ptr) :: bs_brdf_factors ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_brdf_factors_f_shapes
  integer(c_int) :: bs_brdf_factors_f_byte_size

  type(c_ptr) :: bs_nstreams_brdf ! scalar
  integer(c_int) :: bs_nstreams_brdf_f_byte_size

  type(c_ptr) :: bs_do_shadow_effect ! scalar
  integer(c_int) :: bs_do_shadow_effect_f_byte_size

  type(c_ptr) :: bs_do_directbounce_only ! scalar
  integer(c_int) :: bs_do_directbounce_only_f_byte_size

  type(c_ptr) :: bs_do_wsabsa_output ! scalar
  integer(c_int) :: bs_do_wsabsa_output_f_byte_size

  type(c_ptr) :: bs_do_wsa_scaling ! scalar
  integer(c_int) :: bs_do_wsa_scaling_f_byte_size

  type(c_ptr) :: bs_do_bsa_scaling ! scalar
  integer(c_int) :: bs_do_bsa_scaling_f_byte_size

  type(c_ptr) :: bs_wsa_value ! scalar
  integer(c_int) :: bs_wsa_value_f_byte_size

  type(c_ptr) :: bs_bsa_value ! scalar
  integer(c_int) :: bs_bsa_value_f_byte_size

  type(c_ptr) :: bs_do_newcmglint ! scalar
  integer(c_int) :: bs_do_newcmglint_f_byte_size

  type(c_ptr) :: bs_salinity ! scalar
  integer(c_int) :: bs_salinity_f_byte_size

  type(c_ptr) :: bs_wavelength ! scalar
  integer(c_int) :: bs_wavelength_f_byte_size

  type(c_ptr) :: bs_windspeed ! scalar
  integer(c_int) :: bs_windspeed_f_byte_size

  type(c_ptr) :: bs_winddir ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: bs_winddir_f_shapes
  integer(c_int) :: bs_winddir_f_byte_size

  type(c_ptr) :: bs_do_glintshadow ! scalar
  integer(c_int) :: bs_do_glintshadow_f_byte_size

  type(c_ptr) :: bs_do_foamoption ! scalar
  integer(c_int) :: bs_do_foamoption_f_byte_size

  type(c_ptr) :: bs_do_facetisotropy ! scalar
  integer(c_int) :: bs_do_facetisotropy_f_byte_size

  type(c_ptr) :: bs_do_glitter_msrcorr ! scalar
  integer(c_int) :: bs_do_glitter_msrcorr_f_byte_size

  type(c_ptr) :: bs_do_glitter_msrcorr_dbonly ! scalar
  integer(c_int) :: bs_do_glitter_msrcorr_dbonly_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_order ! scalar
  integer(c_int) :: bs_glitter_msrcorr_order_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_nmuquad ! scalar
  integer(c_int) :: bs_glitter_msrcorr_nmuquad_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_nphiquad ! scalar
  integer(c_int) :: bs_glitter_msrcorr_nphiquad_f_byte_size

  
end type brdf_sup_inputs_c

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
type, bind(c) :: brdf_sup_outputs_c
  type(c_ptr) :: bs_dbounce_brdfunc ! dimension(MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(3) :: bs_dbounce_brdfunc_f_shapes
  integer(c_int) :: bs_dbounce_brdfunc_f_byte_size

  type(c_ptr) :: bs_brdf_f_0 ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: bs_brdf_f_0_f_shapes
  integer(c_int) :: bs_brdf_f_0_f_byte_size

  type(c_ptr) :: bs_brdf_f ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS)
  integer(c_int), dimension(3) :: bs_brdf_f_f_shapes
  integer(c_int) :: bs_brdf_f_f_byte_size

  type(c_ptr) :: bs_user_brdf_f_0 ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: bs_user_brdf_f_0_f_shapes
  integer(c_int) :: bs_user_brdf_f_0_f_byte_size

  type(c_ptr) :: bs_user_brdf_f ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS)
  integer(c_int), dimension(3) :: bs_user_brdf_f_f_shapes
  integer(c_int) :: bs_user_brdf_f_f_byte_size

  type(c_ptr) :: bs_emissivity ! dimension(MAXSTREAMS)
  integer(c_int), dimension(1) :: bs_emissivity_f_shapes
  integer(c_int) :: bs_emissivity_f_byte_size

  type(c_ptr) :: bs_user_emissivity ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: bs_user_emissivity_f_shapes
  integer(c_int) :: bs_user_emissivity_f_byte_size

  type(c_ptr) :: bs_wsa_calculated ! scalar
  integer(c_int) :: bs_wsa_calculated_f_byte_size

  type(c_ptr) :: bs_wsa_kernels ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_wsa_kernels_f_shapes
  integer(c_int) :: bs_wsa_kernels_f_byte_size

  type(c_ptr) :: bs_bsa_calculated ! scalar
  integer(c_int) :: bs_bsa_calculated_f_byte_size

  type(c_ptr) :: bs_bsa_kernels ! dimension(MAX_BRDF_KERNELS)
  integer(c_int), dimension(1) :: bs_bsa_kernels_f_shapes
  integer(c_int) :: bs_bsa_kernels_f_byte_size

  
end type brdf_sup_outputs_c

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
type, bind(c) :: brdf_input_exception_handling_c
  type(c_ptr) :: bs_status_inputread ! scalar
  integer(c_int) :: bs_status_inputread_f_byte_size

  type(c_ptr) :: bs_ninputmessages ! scalar
  integer(c_int) :: bs_ninputmessages_f_byte_size

  
  integer(c_int), dimension(1) :: bs_inputmessages_f_shapes
  integer(c_int) :: bs_inputmessages_f_len

  
  integer(c_int), dimension(1) :: bs_inputactions_f_shapes
  integer(c_int) :: bs_inputactions_f_len

  
end type brdf_input_exception_handling_c

! Links to type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
type, bind(c) :: brdf_output_exception_handling_c
  type(c_ptr) :: bs_status_output ! scalar
  integer(c_int) :: bs_status_output_f_byte_size

  type(c_ptr) :: bs_noutputmessages ! scalar
  integer(c_int) :: bs_noutputmessages_f_byte_size

  
  integer(c_int), dimension(1) :: bs_outputmessages_f_shapes
  integer(c_int) :: bs_outputmessages_f_len

  
end type brdf_output_exception_handling_c

! Links to type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.F90"
type, bind(c) :: sleave_sup_inputs_c
  type(c_ptr) :: sl_do_sleaving ! scalar
  integer(c_int) :: sl_do_sleaving_f_byte_size

  type(c_ptr) :: sl_do_isotropic ! scalar
  integer(c_int) :: sl_do_isotropic_f_byte_size

  type(c_ptr) :: sl_do_roughsurface ! scalar
  integer(c_int) :: sl_do_roughsurface_f_byte_size

  type(c_ptr) :: sl_do_exact ! scalar
  integer(c_int) :: sl_do_exact_f_byte_size

  type(c_ptr) :: sl_do_exactonly ! scalar
  integer(c_int) :: sl_do_exactonly_f_byte_size

  type(c_ptr) :: sl_do_fluorescence ! scalar
  integer(c_int) :: sl_do_fluorescence_f_byte_size

  type(c_ptr) :: sl_do_solar_sources ! scalar
  integer(c_int) :: sl_do_solar_sources_f_byte_size

  
  integer(c_int) :: sl_sleave_datapath_f_len

  type(c_ptr) :: sl_do_user_streams ! scalar
  integer(c_int) :: sl_do_user_streams_f_byte_size

  type(c_ptr) :: sl_do_user_obsgeoms ! scalar
  integer(c_int) :: sl_do_user_obsgeoms_f_byte_size

  type(c_ptr) :: sl_do_doublet_geometry ! scalar
  integer(c_int) :: sl_do_doublet_geometry_f_byte_size

  type(c_ptr) :: sl_nstreams ! scalar
  integer(c_int) :: sl_nstreams_f_byte_size

  type(c_ptr) :: sl_nbeams ! scalar
  integer(c_int) :: sl_nbeams_f_byte_size

  type(c_ptr) :: sl_beam_szas ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: sl_beam_szas_f_shapes
  integer(c_int) :: sl_beam_szas_f_byte_size

  type(c_ptr) :: sl_n_user_relazms ! scalar
  integer(c_int) :: sl_n_user_relazms_f_byte_size

  type(c_ptr) :: sl_user_relazms ! dimension(MAX_USER_RELAZMS)
  integer(c_int), dimension(1) :: sl_user_relazms_f_shapes
  integer(c_int) :: sl_user_relazms_f_byte_size

  type(c_ptr) :: sl_n_user_streams ! scalar
  integer(c_int) :: sl_n_user_streams_f_byte_size

  type(c_ptr) :: sl_user_angles_input ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: sl_user_angles_input_f_shapes
  integer(c_int) :: sl_user_angles_input_f_byte_size

  type(c_ptr) :: sl_n_user_obsgeoms ! scalar
  integer(c_int) :: sl_n_user_obsgeoms_f_byte_size

  type(c_ptr) :: sl_user_obsgeoms ! dimension(MAX_USER_OBSGEOMS, 3)
  integer(c_int), dimension(2) :: sl_user_obsgeoms_f_shapes
  integer(c_int) :: sl_user_obsgeoms_f_byte_size

  type(c_ptr) :: sl_n_user_doublets ! scalar
  integer(c_int) :: sl_n_user_doublets_f_byte_size

  type(c_ptr) :: sl_user_doublets ! dimension(MAX_USER_STREAMS, 2)
  integer(c_int), dimension(2) :: sl_user_doublets_f_shapes
  integer(c_int) :: sl_user_doublets_f_byte_size

  type(c_ptr) :: sl_salinity ! scalar
  integer(c_int) :: sl_salinity_f_byte_size

  type(c_ptr) :: sl_chlorconc ! scalar
  integer(c_int) :: sl_chlorconc_f_byte_size

  type(c_ptr) :: sl_wavelength ! scalar
  integer(c_int) :: sl_wavelength_f_byte_size

  type(c_ptr) :: sl_azimuthdep ! scalar
  integer(c_int) :: sl_azimuthdep_f_byte_size

  type(c_ptr) :: sl_do_fourier_output ! scalar
  integer(c_int) :: sl_do_fourier_output_f_byte_size

  type(c_ptr) :: sl_windspeed ! scalar
  integer(c_int) :: sl_windspeed_f_byte_size

  type(c_ptr) :: sl_winddir ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: sl_winddir_f_shapes
  integer(c_int) :: sl_winddir_f_byte_size

  type(c_ptr) :: sl_do_glintshadow ! scalar
  integer(c_int) :: sl_do_glintshadow_f_byte_size

  type(c_ptr) :: sl_do_foamoption ! scalar
  integer(c_int) :: sl_do_foamoption_f_byte_size

  type(c_ptr) :: sl_do_facetisotropy ! scalar
  integer(c_int) :: sl_do_facetisotropy_f_byte_size

  type(c_ptr) :: sl_fl_wavelength ! scalar
  integer(c_int) :: sl_fl_wavelength_f_byte_size

  type(c_ptr) :: sl_fl_latitude ! scalar
  integer(c_int) :: sl_fl_latitude_f_byte_size

  type(c_ptr) :: sl_fl_longitude ! scalar
  integer(c_int) :: sl_fl_longitude_f_byte_size

  type(c_ptr) :: sl_fl_epoch ! dimension(6)
  integer(c_int), dimension(1) :: sl_fl_epoch_f_shapes
  integer(c_int) :: sl_fl_epoch_f_byte_size

  type(c_ptr) :: sl_fl_amplitude755 ! scalar
  integer(c_int) :: sl_fl_amplitude755_f_byte_size

  type(c_ptr) :: sl_fl_do_datagaussian ! scalar
  integer(c_int) :: sl_fl_do_datagaussian_f_byte_size

  type(c_ptr) :: sl_fl_inputgaussians ! dimension(3, 2)
  integer(c_int), dimension(2) :: sl_fl_inputgaussians_f_shapes
  integer(c_int) :: sl_fl_inputgaussians_f_byte_size

  
end type sleave_sup_inputs_c

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_lincontrol_c
  type(c_ptr) :: ts_layer_vary_flag ! dimension(MAXLAYERS)
  integer(c_int), dimension(1) :: ts_layer_vary_flag_f_shapes
  integer(c_int) :: ts_layer_vary_flag_f_byte_size

  type(c_ptr) :: ts_layer_vary_number ! dimension(MAXLAYERS)
  integer(c_int), dimension(1) :: ts_layer_vary_number_f_shapes
  integer(c_int) :: ts_layer_vary_number_f_byte_size

  type(c_ptr) :: ts_n_totalcolumn_wfs ! scalar
  integer(c_int) :: ts_n_totalcolumn_wfs_f_byte_size

  type(c_ptr) :: ts_n_surface_wfs ! scalar
  integer(c_int) :: ts_n_surface_wfs_f_byte_size

  type(c_ptr) :: ts_n_sleave_wfs ! scalar
  integer(c_int) :: ts_n_sleave_wfs_f_byte_size

  
  integer(c_int), dimension(1) :: ts_columnwf_names_f_shapes
  integer(c_int) :: ts_columnwf_names_f_len

  
  integer(c_int), dimension(1) :: ts_profilewf_names_f_shapes
  integer(c_int) :: ts_profilewf_names_f_len

  
end type lidort_fixed_lincontrol_c

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_linoptical_c
  type(c_ptr) :: ts_l_deltau_vert_input ! dimension(MAX_ATMOSWFS, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_l_deltau_vert_input_f_shapes
  integer(c_int) :: ts_l_deltau_vert_input_f_byte_size

  type(c_ptr) :: ts_l_omega_total_input ! dimension(MAX_ATMOSWFS, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_l_omega_total_input_f_shapes
  integer(c_int) :: ts_l_omega_total_input_f_byte_size

  type(c_ptr) :: ts_l_phasmoms_total_input ! dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS)
  integer(c_int), dimension(3) :: ts_l_phasmoms_total_input_f_shapes
  integer(c_int) :: ts_l_phasmoms_total_input_f_byte_size

  type(c_ptr) :: ts_l_phasfunc_input_up ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES)
  integer(c_int), dimension(3) :: ts_l_phasfunc_input_up_f_shapes
  integer(c_int) :: ts_l_phasfunc_input_up_f_byte_size

  type(c_ptr) :: ts_l_phasfunc_input_dn ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES)
  integer(c_int), dimension(3) :: ts_l_phasfunc_input_dn_f_shapes
  integer(c_int) :: ts_l_phasfunc_input_dn_f_byte_size

  
end type lidort_fixed_linoptical_c

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_lininputs_c
  type(c_ptr) :: cont ! scalar
  integer(c_int) :: cont_f_byte_size

  type(c_ptr) :: optical ! scalar
  integer(c_int) :: optical_f_byte_size

  
end type lidort_fixed_lininputs_c

! Links to type: "lidort_modified_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_modified_lincontrol_c
  type(c_ptr) :: ts_do_column_linearization ! scalar
  integer(c_int) :: ts_do_column_linearization_f_byte_size

  type(c_ptr) :: ts_do_profile_linearization ! scalar
  integer(c_int) :: ts_do_profile_linearization_f_byte_size

  type(c_ptr) :: ts_do_atmos_linearization ! scalar
  integer(c_int) :: ts_do_atmos_linearization_f_byte_size

  type(c_ptr) :: ts_do_surface_linearization ! scalar
  integer(c_int) :: ts_do_surface_linearization_f_byte_size

  type(c_ptr) :: ts_do_linearization ! scalar
  integer(c_int) :: ts_do_linearization_f_byte_size

  type(c_ptr) :: ts_do_simulation_only ! scalar
  integer(c_int) :: ts_do_simulation_only_f_byte_size

  type(c_ptr) :: ts_do_atmos_lbbf ! scalar
  integer(c_int) :: ts_do_atmos_lbbf_f_byte_size

  type(c_ptr) :: ts_do_surface_lbbf ! scalar
  integer(c_int) :: ts_do_surface_lbbf_f_byte_size

  type(c_ptr) :: ts_do_sleave_wfs ! scalar
  integer(c_int) :: ts_do_sleave_wfs_f_byte_size

  
end type lidort_modified_lincontrol_c

! Links to type: "lidort_modified_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_modified_lininputs_c
  type(c_ptr) :: mcont ! scalar
  integer(c_int) :: mcont_f_byte_size

  
end type lidort_modified_lininputs_c

! Links to type: "lidort_linatmos" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linatmos_c
  type(c_ptr) :: ts_columnwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_columnwf_f_shapes
  integer(c_int) :: ts_columnwf_f_byte_size

  type(c_ptr) :: ts_meani_diffuse_colwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_meani_diffuse_colwf_f_shapes
  integer(c_int) :: ts_meani_diffuse_colwf_f_byte_size

  type(c_ptr) :: ts_flux_diffuse_colwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_flux_diffuse_colwf_f_shapes
  integer(c_int) :: ts_flux_diffuse_colwf_f_byte_size

  type(c_ptr) :: ts_dnmeani_direct_colwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_dnmeani_direct_colwf_f_shapes
  integer(c_int) :: ts_dnmeani_direct_colwf_f_byte_size

  type(c_ptr) :: ts_dnflux_direct_colwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_dnflux_direct_colwf_f_shapes
  integer(c_int) :: ts_dnflux_direct_colwf_f_byte_size

  type(c_ptr) :: ts_profilewf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(5) :: ts_profilewf_f_shapes
  integer(c_int) :: ts_profilewf_f_byte_size

  type(c_ptr) :: ts_meani_diffuse_profwf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(5) :: ts_meani_diffuse_profwf_f_shapes
  integer(c_int) :: ts_meani_diffuse_profwf_f_byte_size

  type(c_ptr) :: ts_flux_diffuse_profwf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(5) :: ts_flux_diffuse_profwf_f_shapes
  integer(c_int) :: ts_flux_diffuse_profwf_f_byte_size

  type(c_ptr) :: ts_dnmeani_direct_profwf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_dnmeani_direct_profwf_f_shapes
  integer(c_int) :: ts_dnmeani_direct_profwf_f_byte_size

  type(c_ptr) :: ts_dnflux_direct_profwf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_dnflux_direct_profwf_f_shapes
  integer(c_int) :: ts_dnflux_direct_profwf_f_byte_size

  type(c_ptr) :: ts_abbwfs_jacobians ! dimension(MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_abbwfs_jacobians_f_shapes
  integer(c_int) :: ts_abbwfs_jacobians_f_byte_size

  type(c_ptr) :: ts_abbwfs_fluxes ! dimension(MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_abbwfs_fluxes_f_shapes
  integer(c_int) :: ts_abbwfs_fluxes_f_byte_size

  type(c_ptr) :: ts_albmed_user_profwf ! dimension(MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_albmed_user_profwf_f_shapes
  integer(c_int) :: ts_albmed_user_profwf_f_byte_size

  type(c_ptr) :: ts_trnmed_user_profwf ! dimension(MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_trnmed_user_profwf_f_shapes
  integer(c_int) :: ts_trnmed_user_profwf_f_byte_size

  type(c_ptr) :: ts_albmed_fluxes_profwf ! dimension(2, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_albmed_fluxes_profwf_f_shapes
  integer(c_int) :: ts_albmed_fluxes_profwf_f_byte_size

  type(c_ptr) :: ts_trnmed_fluxes_profwf ! dimension(2, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_trnmed_fluxes_profwf_f_shapes
  integer(c_int) :: ts_trnmed_fluxes_profwf_f_byte_size

  type(c_ptr) :: ts_transbeam_profwf ! dimension(MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_transbeam_profwf_f_shapes
  integer(c_int) :: ts_transbeam_profwf_f_byte_size

  type(c_ptr) :: ts_albmed_user_colwf ! dimension(MAX_USER_STREAMS, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_albmed_user_colwf_f_shapes
  integer(c_int) :: ts_albmed_user_colwf_f_byte_size

  type(c_ptr) :: ts_trnmed_user_colwf ! dimension(MAX_USER_STREAMS, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_trnmed_user_colwf_f_shapes
  integer(c_int) :: ts_trnmed_user_colwf_f_byte_size

  type(c_ptr) :: ts_albmed_fluxes_colwf ! dimension(2, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_albmed_fluxes_colwf_f_shapes
  integer(c_int) :: ts_albmed_fluxes_colwf_f_byte_size

  type(c_ptr) :: ts_trnmed_fluxes_colwf ! dimension(2, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_trnmed_fluxes_colwf_f_shapes
  integer(c_int) :: ts_trnmed_fluxes_colwf_f_byte_size

  type(c_ptr) :: ts_transbeam_colwf ! dimension(MAXBEAMS, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_transbeam_colwf_f_shapes
  integer(c_int) :: ts_transbeam_colwf_f_byte_size

  type(c_ptr) :: ts_planetary_transterm_profwf ! dimension(MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(3) :: ts_planetary_transterm_profwf_f_shapes
  integer(c_int) :: ts_planetary_transterm_profwf_f_byte_size

  type(c_ptr) :: ts_planetary_sbterm_profwf ! dimension(MAXLAYERS, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_planetary_sbterm_profwf_f_shapes
  integer(c_int) :: ts_planetary_sbterm_profwf_f_byte_size

  type(c_ptr) :: ts_planetary_transterm_colwf ! dimension(MAX_GEOMETRIES, MAX_ATMOSWFS)
  integer(c_int), dimension(2) :: ts_planetary_transterm_colwf_f_shapes
  integer(c_int) :: ts_planetary_transterm_colwf_f_byte_size

  type(c_ptr) :: ts_planetary_sbterm_colwf ! dimension(MAX_ATMOSWFS)
  integer(c_int), dimension(1) :: ts_planetary_sbterm_colwf_f_shapes
  integer(c_int) :: ts_planetary_sbterm_colwf_f_byte_size

  type(c_ptr) :: ts_lc_lostrans ! dimension(MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(3) :: ts_lc_lostrans_f_shapes
  integer(c_int) :: ts_lc_lostrans_f_byte_size

  type(c_ptr) :: ts_lc_layer_mssts ! dimension(MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(3) :: ts_lc_layer_mssts_f_shapes
  integer(c_int) :: ts_lc_layer_mssts_f_byte_size

  type(c_ptr) :: ts_lc_surf_mssts ! dimension(MAX_ATMOSWFS, MAXBEAMS)
  integer(c_int), dimension(2) :: ts_lc_surf_mssts_f_shapes
  integer(c_int) :: ts_lc_surf_mssts_f_byte_size

  type(c_ptr) :: ts_lp_lostrans ! dimension(MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(3) :: ts_lp_lostrans_f_shapes
  integer(c_int) :: ts_lp_lostrans_f_byte_size

  type(c_ptr) :: ts_lp_layer_mssts ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(4) :: ts_lp_layer_mssts_f_shapes
  integer(c_int) :: ts_lp_layer_mssts_f_byte_size

  type(c_ptr) :: ts_lp_surf_mssts ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_lp_surf_mssts_f_shapes
  integer(c_int) :: ts_lp_surf_mssts_f_byte_size

  
end type lidort_linatmos_c

! Links to type: "lidort_linsurf" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linsurf_c
  type(c_ptr) :: ts_surfacewf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_surfacewf_f_shapes
  integer(c_int) :: ts_surfacewf_f_byte_size

  type(c_ptr) :: ts_meani_diffuse_surfwf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_meani_diffuse_surfwf_f_shapes
  integer(c_int) :: ts_meani_diffuse_surfwf_f_byte_size

  type(c_ptr) :: ts_flux_diffuse_surfwf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_flux_diffuse_surfwf_f_shapes
  integer(c_int) :: ts_flux_diffuse_surfwf_f_byte_size

  type(c_ptr) :: ts_sbbwfs_jacobians ! dimension(MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_sbbwfs_jacobians_f_shapes
  integer(c_int) :: ts_sbbwfs_jacobians_f_byte_size

  type(c_ptr) :: ts_sbbwfs_fluxes ! dimension(MAX_USER_LEVELS, 2, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_sbbwfs_fluxes_f_shapes
  integer(c_int) :: ts_sbbwfs_fluxes_f_byte_size

  type(c_ptr) :: ts_ls_layer_mssts ! dimension(MAX_SURFACEWFS, MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(3) :: ts_ls_layer_mssts_f_shapes
  integer(c_int) :: ts_ls_layer_mssts_f_byte_size

  type(c_ptr) :: ts_ls_surf_mssts ! dimension(MAX_SURFACEWFS, MAXBEAMS)
  integer(c_int), dimension(2) :: ts_ls_surf_mssts_f_shapes
  integer(c_int) :: ts_ls_surf_mssts_f_byte_size

  
end type lidort_linsurf_c

! Links to type: "lidort_linoutputs" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linoutputs_c
  type(c_ptr) :: atmos ! scalar
  integer(c_int) :: atmos_f_byte_size

  type(c_ptr) :: surf ! scalar
  integer(c_int) :: surf_f_byte_size

  
end type lidort_linoutputs_c

! Links to type: "lidort_linsup_brdf" from module: "lidort_lin_sup_brdf_def_m" in file: "lidort_lin_sup_brdf_def.F90"
type, bind(c) :: lidort_linsup_brdf_c
  type(c_ptr) :: ts_ls_exactdb_brdfunc ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_ls_exactdb_brdfunc_f_shapes
  integer(c_int) :: ts_ls_exactdb_brdfunc_f_byte_size

  type(c_ptr) :: ts_ls_brdf_f_0 ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_ls_brdf_f_0_f_shapes
  integer(c_int) :: ts_ls_brdf_f_0_f_byte_size

  type(c_ptr) :: ts_ls_brdf_f ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS)
  integer(c_int), dimension(4) :: ts_ls_brdf_f_f_shapes
  integer(c_int) :: ts_ls_brdf_f_f_byte_size

  type(c_ptr) :: ts_ls_user_brdf_f_0 ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_ls_user_brdf_f_0_f_shapes
  integer(c_int) :: ts_ls_user_brdf_f_0_f_byte_size

  type(c_ptr) :: ts_ls_user_brdf_f ! dimension(MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS)
  integer(c_int), dimension(4) :: ts_ls_user_brdf_f_f_shapes
  integer(c_int) :: ts_ls_user_brdf_f_f_byte_size

  type(c_ptr) :: ts_ls_emissivity ! dimension(MAX_SURFACEWFS, MAXSTREAMS)
  integer(c_int), dimension(2) :: ts_ls_emissivity_f_shapes
  integer(c_int) :: ts_ls_emissivity_f_byte_size

  type(c_ptr) :: ts_ls_user_emissivity ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS)
  integer(c_int), dimension(2) :: ts_ls_user_emissivity_f_shapes
  integer(c_int) :: ts_ls_user_emissivity_f_byte_size

  
end type lidort_linsup_brdf_c

! Links to type: "lidort_linsup_sleave" from module: "lidort_lin_sup_sleave_def_m" in file: "lidort_lin_sup_sleave_def.F90"
type, bind(c) :: lidort_linsup_sleave_c
  type(c_ptr) :: ts_lssl_slterm_isotropic ! dimension(MAX_SLEAVEWFS, MAXBEAMS)
  integer(c_int), dimension(2) :: ts_lssl_slterm_isotropic_f_shapes
  integer(c_int) :: ts_lssl_slterm_isotropic_f_byte_size

  type(c_ptr) :: ts_lssl_slterm_userangles ! dimension(MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_lssl_slterm_userangles_f_shapes
  integer(c_int) :: ts_lssl_slterm_userangles_f_byte_size

  type(c_ptr) :: ts_lssl_slterm_f_0 ! dimension(MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_lssl_slterm_f_0_f_shapes
  integer(c_int) :: ts_lssl_slterm_f_0_f_byte_size

  type(c_ptr) :: ts_lssl_user_slterm_f_0 ! dimension(MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(4) :: ts_lssl_user_slterm_f_0_f_shapes
  integer(c_int) :: ts_lssl_user_slterm_f_0_f_byte_size

  
end type lidort_linsup_sleave_c

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
type, bind(c) :: lidort_linsup_ss_atmos_c
  type(c_ptr) :: ts_columnwf_ss ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(4) :: ts_columnwf_ss_f_shapes
  integer(c_int) :: ts_columnwf_ss_f_byte_size

  type(c_ptr) :: ts_columnwf_db ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(3) :: ts_columnwf_db_f_shapes
  integer(c_int) :: ts_columnwf_db_f_byte_size

  type(c_ptr) :: ts_profilewf_ss ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(5) :: ts_profilewf_ss_f_shapes
  integer(c_int) :: ts_profilewf_ss_f_byte_size

  type(c_ptr) :: ts_profilewf_db ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(4) :: ts_profilewf_db_f_shapes
  integer(c_int) :: ts_profilewf_db_f_byte_size

  
end type lidort_linsup_ss_atmos_c

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
type, bind(c) :: lidort_linsup_ss_surf_c
  type(c_ptr) :: ts_surfacewf_db ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(3) :: ts_surfacewf_db_f_shapes
  integer(c_int) :: ts_surfacewf_db_f_byte_size

  
end type lidort_linsup_ss_surf_c

! Links to type: "lidort_linsup_ss" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
type, bind(c) :: lidort_linsup_ss_c
  type(c_ptr) :: atmos ! scalar
  integer(c_int) :: atmos_f_byte_size

  type(c_ptr) :: surf ! scalar
  integer(c_int) :: surf_f_byte_size

  
end type lidort_linsup_ss_c

! Links to type: "lidort_linsup_inout" from module: "lidort_lin_sup_inout_def_m" in file: "lidort_lin_sup_def.F90"
type, bind(c) :: lidort_linsup_inout_c
  type(c_ptr) :: brdf ! scalar
  integer(c_int) :: brdf_f_byte_size

  type(c_ptr) :: ss ! scalar
  integer(c_int) :: ss_f_byte_size

  type(c_ptr) :: sleave ! scalar
  integer(c_int) :: sleave_f_byte_size

  
end type lidort_linsup_inout_c

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_main_outputs_c
  type(c_ptr) :: ts_intensity ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_intensity_f_shapes
  integer(c_int) :: ts_intensity_f_byte_size

  type(c_ptr) :: ts_meani_diffuse ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_meani_diffuse_f_shapes
  integer(c_int) :: ts_meani_diffuse_f_byte_size

  type(c_ptr) :: ts_flux_diffuse ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_flux_diffuse_f_shapes
  integer(c_int) :: ts_flux_diffuse_f_byte_size

  type(c_ptr) :: ts_dnmeani_direct ! dimension(MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(2) :: ts_dnmeani_direct_f_shapes
  integer(c_int) :: ts_dnmeani_direct_f_byte_size

  type(c_ptr) :: ts_dnflux_direct ! dimension(MAX_USER_LEVELS, MAXBEAMS)
  integer(c_int), dimension(2) :: ts_dnflux_direct_f_shapes
  integer(c_int) :: ts_dnflux_direct_f_byte_size

  type(c_ptr) :: ts_albmed_user ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_albmed_user_f_shapes
  integer(c_int) :: ts_albmed_user_f_byte_size

  type(c_ptr) :: ts_trnmed_user ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_trnmed_user_f_shapes
  integer(c_int) :: ts_trnmed_user_f_byte_size

  type(c_ptr) :: ts_albmed_fluxes ! dimension(2)
  integer(c_int), dimension(1) :: ts_albmed_fluxes_f_shapes
  integer(c_int) :: ts_albmed_fluxes_f_byte_size

  type(c_ptr) :: ts_trnmed_fluxes ! dimension(2)
  integer(c_int), dimension(1) :: ts_trnmed_fluxes_f_shapes
  integer(c_int) :: ts_trnmed_fluxes_f_byte_size

  type(c_ptr) :: ts_planetary_transterm ! dimension(MAX_GEOMETRIES)
  integer(c_int), dimension(1) :: ts_planetary_transterm_f_shapes
  integer(c_int) :: ts_planetary_transterm_f_byte_size

  type(c_ptr) :: ts_planetary_sbterm ! scalar
  integer(c_int) :: ts_planetary_sbterm_f_byte_size

  type(c_ptr) :: ts_pathgeoms ! dimension(2, 0:MAXLAYERS)
  integer(c_int), dimension(2) :: ts_pathgeoms_f_shapes
  integer(c_int) :: ts_pathgeoms_f_byte_size

  type(c_ptr) :: ts_lostrans ! dimension(MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_lostrans_f_shapes
  integer(c_int) :: ts_lostrans_f_byte_size

  type(c_ptr) :: ts_layer_mssts ! dimension(MAXBEAMS, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_layer_mssts_f_shapes
  integer(c_int) :: ts_layer_mssts_f_byte_size

  type(c_ptr) :: ts_surf_mssts ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_surf_mssts_f_shapes
  integer(c_int) :: ts_surf_mssts_f_byte_size

  type(c_ptr) :: ts_contribs ! dimension(MAX_GEOMETRIES, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_contribs_f_shapes
  integer(c_int) :: ts_contribs_f_byte_size

  type(c_ptr) :: ts_fourier_saved ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_fourier_saved_f_shapes
  integer(c_int) :: ts_fourier_saved_f_byte_size

  type(c_ptr) :: ts_n_geometries ! scalar
  integer(c_int) :: ts_n_geometries_f_byte_size

  type(c_ptr) :: ts_solarbeam_boatrans ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_solarbeam_boatrans_f_shapes
  integer(c_int) :: ts_solarbeam_boatrans_f_byte_size

  type(c_ptr) :: ts_spheralb ! scalar
  integer(c_int) :: ts_spheralb_f_byte_size

  type(c_ptr) :: ts_trans1_user ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_trans1_user_f_shapes
  integer(c_int) :: ts_trans1_user_f_byte_size

  type(c_ptr) :: ts_trans1_beam ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_trans1_beam_f_shapes
  integer(c_int) :: ts_trans1_beam_f_byte_size

  
end type lidort_main_outputs_c

! Links to type: "lidort_wladjusted_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_wladjusted_outputs_c
  type(c_ptr) :: ts_wladjusted_isotropic ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_wladjusted_isotropic_f_shapes
  integer(c_int) :: ts_wladjusted_isotropic_f_byte_size

  type(c_ptr) :: ts_wladjusted_direct ! dimension(MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_wladjusted_direct_f_shapes
  integer(c_int) :: ts_wladjusted_direct_f_byte_size

  type(c_ptr) :: ts_wladjusted_f_ords_0 ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_wladjusted_f_ords_0_f_shapes
  integer(c_int) :: ts_wladjusted_f_ords_0_f_byte_size

  type(c_ptr) :: ts_wladjusted_f_user_0 ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_wladjusted_f_user_0_f_shapes
  integer(c_int) :: ts_wladjusted_f_user_0_f_byte_size

  
end type lidort_wladjusted_outputs_c

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_exception_handling_c
  type(c_ptr) :: ts_status_inputcheck ! scalar
  integer(c_int) :: ts_status_inputcheck_f_byte_size

  type(c_ptr) :: ts_ncheckmessages ! scalar
  integer(c_int) :: ts_ncheckmessages_f_byte_size

  
  integer(c_int), dimension(1) :: ts_checkmessages_f_shapes
  integer(c_int) :: ts_checkmessages_f_len

  
  integer(c_int), dimension(1) :: ts_actions_f_shapes
  integer(c_int) :: ts_actions_f_len

  type(c_ptr) :: ts_status_calculation ! scalar
  integer(c_int) :: ts_status_calculation_f_byte_size

  
  integer(c_int) :: ts_message_f_len

  
  integer(c_int) :: ts_trace_1_f_len

  
  integer(c_int) :: ts_trace_2_f_len

  
  integer(c_int) :: ts_trace_3_f_len

  
end type lidort_exception_handling_c

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_input_exception_handling_c
  type(c_ptr) :: ts_status_inputread ! scalar
  integer(c_int) :: ts_status_inputread_f_byte_size

  type(c_ptr) :: ts_ninputmessages ! scalar
  integer(c_int) :: ts_ninputmessages_f_byte_size

  
  integer(c_int), dimension(1) :: ts_inputmessages_f_shapes
  integer(c_int) :: ts_inputmessages_f_len

  
  integer(c_int), dimension(1) :: ts_inputactions_f_shapes
  integer(c_int) :: ts_inputactions_f_len

  
end type lidort_input_exception_handling_c

! Links to type: "lidort_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_outputs_c
  type(c_ptr) :: main ! scalar
  integer(c_int) :: main_f_byte_size

  type(c_ptr) :: wlout ! scalar
  integer(c_int) :: wlout_f_byte_size

  type(c_ptr) :: status ! scalar
  integer(c_int) :: status_f_byte_size

  
end type lidort_outputs_c

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def_m" in file: "lidort_sup_brdf_def.F90"
type, bind(c) :: lidort_sup_brdf_c
  type(c_ptr) :: ts_exactdb_brdfunc ! dimension(MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_exactdb_brdfunc_f_shapes
  integer(c_int) :: ts_exactdb_brdfunc_f_byte_size

  type(c_ptr) :: ts_brdf_f_0 ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_brdf_f_0_f_shapes
  integer(c_int) :: ts_brdf_f_0_f_byte_size

  type(c_ptr) :: ts_brdf_f ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS)
  integer(c_int), dimension(3) :: ts_brdf_f_f_shapes
  integer(c_int) :: ts_brdf_f_f_byte_size

  type(c_ptr) :: ts_user_brdf_f_0 ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_user_brdf_f_0_f_shapes
  integer(c_int) :: ts_user_brdf_f_0_f_byte_size

  type(c_ptr) :: ts_user_brdf_f ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS)
  integer(c_int), dimension(3) :: ts_user_brdf_f_f_shapes
  integer(c_int) :: ts_user_brdf_f_f_byte_size

  type(c_ptr) :: ts_emissivity ! dimension(MAXSTREAMS)
  integer(c_int), dimension(1) :: ts_emissivity_f_shapes
  integer(c_int) :: ts_emissivity_f_byte_size

  type(c_ptr) :: ts_user_emissivity ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_user_emissivity_f_shapes
  integer(c_int) :: ts_user_emissivity_f_byte_size

  
end type lidort_sup_brdf_c

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def_m" in file: "lidort_sup_sleave_def.F90"
type, bind(c) :: lidort_sup_sleave_c
  type(c_ptr) :: ts_slterm_isotropic ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_slterm_isotropic_f_shapes
  integer(c_int) :: ts_slterm_isotropic_f_byte_size

  type(c_ptr) :: ts_slterm_userangles ! dimension(MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_slterm_userangles_f_shapes
  integer(c_int) :: ts_slterm_userangles_f_byte_size

  type(c_ptr) :: ts_slterm_f_0 ! dimension(0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_slterm_f_0_f_shapes
  integer(c_int) :: ts_slterm_f_0_f_byte_size

  type(c_ptr) :: ts_user_slterm_f_0 ! dimension(0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS)
  integer(c_int), dimension(3) :: ts_user_slterm_f_0_f_shapes
  integer(c_int) :: ts_user_slterm_f_0_f_byte_size

  
end type lidort_sup_sleave_c

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def_m" in file: "lidort_sup_ss_def.F90"
type, bind(c) :: lidort_sup_ss_c
  type(c_ptr) :: ts_intensity_ss ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_intensity_ss_f_shapes
  integer(c_int) :: ts_intensity_ss_f_byte_size

  type(c_ptr) :: ts_intensity_db ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(2) :: ts_intensity_db_f_shapes
  integer(c_int) :: ts_intensity_db_f_byte_size

  type(c_ptr) :: ts_contribs_ss ! dimension(MAX_GEOMETRIES, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_contribs_ss_f_shapes
  integer(c_int) :: ts_contribs_ss_f_byte_size

  
end type lidort_sup_ss_c

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def_m" in file: "lidort_sup_def.F90"
type, bind(c) :: lidort_sup_inout_c
  type(c_ptr) :: brdf ! scalar
  integer(c_int) :: brdf_f_byte_size

  type(c_ptr) :: ss ! scalar
  integer(c_int) :: ss_f_byte_size

  type(c_ptr) :: sleave ! scalar
  integer(c_int) :: sleave_f_byte_size

  
end type lidort_sup_inout_c

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_boolean_c
  type(c_ptr) :: ts_do_fullrad_mode ! scalar
  integer(c_int) :: ts_do_fullrad_mode_f_byte_size

  type(c_ptr) :: ts_do_thermal_emission ! scalar
  integer(c_int) :: ts_do_thermal_emission_f_byte_size

  type(c_ptr) :: ts_do_surface_emission ! scalar
  integer(c_int) :: ts_do_surface_emission_f_byte_size

  type(c_ptr) :: ts_do_plane_parallel ! scalar
  integer(c_int) :: ts_do_plane_parallel_f_byte_size

  type(c_ptr) :: ts_do_brdf_surface ! scalar
  integer(c_int) :: ts_do_brdf_surface_f_byte_size

  type(c_ptr) :: ts_do_upwelling ! scalar
  integer(c_int) :: ts_do_upwelling_f_byte_size

  type(c_ptr) :: ts_do_dnwelling ! scalar
  integer(c_int) :: ts_do_dnwelling_f_byte_size

  type(c_ptr) :: ts_do_toa_contribs ! scalar
  integer(c_int) :: ts_do_toa_contribs_f_byte_size

  type(c_ptr) :: ts_do_surface_leaving ! scalar
  integer(c_int) :: ts_do_surface_leaving_f_byte_size

  type(c_ptr) :: ts_do_sl_isotropic ! scalar
  integer(c_int) :: ts_do_sl_isotropic_f_byte_size

  type(c_ptr) :: ts_do_water_leaving ! scalar
  integer(c_int) :: ts_do_water_leaving_f_byte_size

  type(c_ptr) :: ts_do_fluorescence ! scalar
  integer(c_int) :: ts_do_fluorescence_f_byte_size

  type(c_ptr) :: ts_do_tf_iteration ! scalar
  integer(c_int) :: ts_do_tf_iteration_f_byte_size

  type(c_ptr) :: ts_do_wladjusted_output ! scalar
  integer(c_int) :: ts_do_wladjusted_output_f_byte_size

  type(c_ptr) :: ts_do_toa_illumination ! scalar
  integer(c_int) :: ts_do_toa_illumination_f_byte_size

  type(c_ptr) :: ts_do_boa_illumination ! scalar
  integer(c_int) :: ts_do_boa_illumination_f_byte_size

  type(c_ptr) :: ts_do_albtrn_media ! dimension(2)
  integer(c_int), dimension(1) :: ts_do_albtrn_media_f_shapes
  integer(c_int) :: ts_do_albtrn_media_f_byte_size

  type(c_ptr) :: ts_do_planetary_problem ! scalar
  integer(c_int) :: ts_do_planetary_problem_f_byte_size

  type(c_ptr) :: ts_do_mssts ! scalar
  integer(c_int) :: ts_do_mssts_f_byte_size

  
end type lidort_fixed_boolean_c

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_control_c
  type(c_ptr) :: ts_taylor_order ! scalar
  integer(c_int) :: ts_taylor_order_f_byte_size

  type(c_ptr) :: ts_nstreams ! scalar
  integer(c_int) :: ts_nstreams_f_byte_size

  type(c_ptr) :: ts_nlayers ! scalar
  integer(c_int) :: ts_nlayers_f_byte_size

  type(c_ptr) :: ts_nfinelayers ! scalar
  integer(c_int) :: ts_nfinelayers_f_byte_size

  type(c_ptr) :: ts_n_thermal_coeffs ! scalar
  integer(c_int) :: ts_n_thermal_coeffs_f_byte_size

  type(c_ptr) :: ts_lidort_accuracy ! scalar
  integer(c_int) :: ts_lidort_accuracy_f_byte_size

  type(c_ptr) :: ts_asymtx_tolerance ! scalar
  integer(c_int) :: ts_asymtx_tolerance_f_byte_size

  type(c_ptr) :: ts_tf_maxiter ! scalar
  integer(c_int) :: ts_tf_maxiter_f_byte_size

  type(c_ptr) :: ts_tf_criterion ! scalar
  integer(c_int) :: ts_tf_criterion_f_byte_size

  type(c_ptr) :: ts_toa_illumination ! scalar
  integer(c_int) :: ts_toa_illumination_f_byte_size

  type(c_ptr) :: ts_boa_illumination ! scalar
  integer(c_int) :: ts_boa_illumination_f_byte_size

  
end type lidort_fixed_control_c

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_sunrays_c
  type(c_ptr) :: ts_flux_factor ! scalar
  integer(c_int) :: ts_flux_factor_f_byte_size

  
end type lidort_fixed_sunrays_c

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_uservalues_c
  type(c_ptr) :: ts_n_user_levels ! scalar
  integer(c_int) :: ts_n_user_levels_f_byte_size

  
end type lidort_fixed_uservalues_c

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_chapman_c
  type(c_ptr) :: ts_height_grid ! dimension(0:MAXLAYERS)
  integer(c_int), dimension(1) :: ts_height_grid_f_shapes
  integer(c_int) :: ts_height_grid_f_byte_size

  type(c_ptr) :: ts_pressure_grid ! dimension(0:MAXLAYERS)
  integer(c_int), dimension(1) :: ts_pressure_grid_f_shapes
  integer(c_int) :: ts_pressure_grid_f_byte_size

  type(c_ptr) :: ts_temperature_grid ! dimension(0:MAXLAYERS)
  integer(c_int), dimension(1) :: ts_temperature_grid_f_shapes
  integer(c_int) :: ts_temperature_grid_f_byte_size

  type(c_ptr) :: ts_finegrid ! dimension(MAXLAYERS)
  integer(c_int), dimension(1) :: ts_finegrid_f_shapes
  integer(c_int) :: ts_finegrid_f_byte_size

  type(c_ptr) :: ts_rfindex_parameter ! scalar
  integer(c_int) :: ts_rfindex_parameter_f_byte_size

  
end type lidort_fixed_chapman_c

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_optical_c
  type(c_ptr) :: ts_deltau_vert_input ! dimension(MAXLAYERS)
  integer(c_int), dimension(1) :: ts_deltau_vert_input_f_shapes
  integer(c_int) :: ts_deltau_vert_input_f_byte_size

  type(c_ptr) :: ts_phasmoms_total_input ! dimension(0:MAXMOMENTS_INPUT, MAXLAYERS)
  integer(c_int), dimension(2) :: ts_phasmoms_total_input_f_shapes
  integer(c_int) :: ts_phasmoms_total_input_f_byte_size

  type(c_ptr) :: ts_phasfunc_input_up ! dimension(MAXLAYERS, MAX_GEOMETRIES)
  integer(c_int), dimension(2) :: ts_phasfunc_input_up_f_shapes
  integer(c_int) :: ts_phasfunc_input_up_f_byte_size

  type(c_ptr) :: ts_phasfunc_input_dn ! dimension(MAXLAYERS, MAX_GEOMETRIES)
  integer(c_int), dimension(2) :: ts_phasfunc_input_dn_f_shapes
  integer(c_int) :: ts_phasfunc_input_dn_f_byte_size

  type(c_ptr) :: ts_lambertian_albedo ! scalar
  integer(c_int) :: ts_lambertian_albedo_f_byte_size

  type(c_ptr) :: ts_thermal_bb_input ! dimension(0:MAXLAYERS)
  integer(c_int), dimension(1) :: ts_thermal_bb_input_f_shapes
  integer(c_int) :: ts_thermal_bb_input_f_byte_size

  type(c_ptr) :: ts_surface_bb_input ! scalar
  integer(c_int) :: ts_surface_bb_input_f_byte_size

  type(c_ptr) :: ts_atmos_wavelength ! scalar
  integer(c_int) :: ts_atmos_wavelength_f_byte_size

  
end type lidort_fixed_optical_c

! Links to type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_write_c
  type(c_ptr) :: ts_do_debug_write ! scalar
  integer(c_int) :: ts_do_debug_write_f_byte_size

  type(c_ptr) :: ts_do_write_input ! scalar
  integer(c_int) :: ts_do_write_input_f_byte_size

  
  integer(c_int) :: ts_input_write_filename_f_len

  type(c_ptr) :: ts_do_write_scenario ! scalar
  integer(c_int) :: ts_do_write_scenario_f_byte_size

  
  integer(c_int) :: ts_scenario_write_filename_f_len

  type(c_ptr) :: ts_do_write_fourier ! scalar
  integer(c_int) :: ts_do_write_fourier_f_byte_size

  
  integer(c_int) :: ts_fourier_write_filename_f_len

  type(c_ptr) :: ts_do_write_results ! scalar
  integer(c_int) :: ts_do_write_results_f_byte_size

  
  integer(c_int) :: ts_results_write_filename_f_len

  
end type lidort_fixed_write_c

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_inputs_c
  type(c_ptr) :: f_bool ! scalar
  integer(c_int) :: bool_f_byte_size

  type(c_ptr) :: cont ! scalar
  integer(c_int) :: cont_f_byte_size

  type(c_ptr) :: sunrays ! scalar
  integer(c_int) :: sunrays_f_byte_size

  type(c_ptr) :: userval ! scalar
  integer(c_int) :: userval_f_byte_size

  type(c_ptr) :: chapman ! scalar
  integer(c_int) :: chapman_f_byte_size

  type(c_ptr) :: optical ! scalar
  integer(c_int) :: optical_f_byte_size

  type(c_ptr) :: write ! scalar
  integer(c_int) :: write_f_byte_size

  
end type lidort_fixed_inputs_c

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_boolean_c
  type(c_ptr) :: ts_do_focorr ! scalar
  integer(c_int) :: ts_do_focorr_f_byte_size

  type(c_ptr) :: ts_do_focorr_external ! scalar
  integer(c_int) :: ts_do_focorr_external_f_byte_size

  type(c_ptr) :: ts_do_focorr_nadir ! scalar
  integer(c_int) :: ts_do_focorr_nadir_f_byte_size

  type(c_ptr) :: ts_do_focorr_outgoing ! scalar
  integer(c_int) :: ts_do_focorr_outgoing_f_byte_size

  type(c_ptr) :: ts_do_sscorr_truncation ! scalar
  integer(c_int) :: ts_do_sscorr_truncation_f_byte_size

  type(c_ptr) :: ts_do_sscorr_usephasfunc ! scalar
  integer(c_int) :: ts_do_sscorr_usephasfunc_f_byte_size

  type(c_ptr) :: ts_do_external_wleave ! scalar
  integer(c_int) :: ts_do_external_wleave_f_byte_size

  type(c_ptr) :: ts_do_double_convtest ! scalar
  integer(c_int) :: ts_do_double_convtest_f_byte_size

  type(c_ptr) :: ts_do_solar_sources ! scalar
  integer(c_int) :: ts_do_solar_sources_f_byte_size

  type(c_ptr) :: ts_do_refractive_geometry ! scalar
  integer(c_int) :: ts_do_refractive_geometry_f_byte_size

  type(c_ptr) :: ts_do_chapman_function ! scalar
  integer(c_int) :: ts_do_chapman_function_f_byte_size

  type(c_ptr) :: ts_do_rayleigh_only ! scalar
  integer(c_int) :: ts_do_rayleigh_only_f_byte_size

  type(c_ptr) :: ts_do_isotropic_only ! scalar
  integer(c_int) :: ts_do_isotropic_only_f_byte_size

  type(c_ptr) :: ts_do_no_azimuth ! scalar
  integer(c_int) :: ts_do_no_azimuth_f_byte_size

  type(c_ptr) :: ts_do_all_fourier ! scalar
  integer(c_int) :: ts_do_all_fourier_f_byte_size

  type(c_ptr) :: ts_do_deltam_scaling ! scalar
  integer(c_int) :: ts_do_deltam_scaling_f_byte_size

  type(c_ptr) :: ts_do_solution_saving ! scalar
  integer(c_int) :: ts_do_solution_saving_f_byte_size

  type(c_ptr) :: ts_do_bvp_telescoping ! scalar
  integer(c_int) :: ts_do_bvp_telescoping_f_byte_size

  type(c_ptr) :: ts_do_user_streams ! scalar
  integer(c_int) :: ts_do_user_streams_f_byte_size

  type(c_ptr) :: ts_do_additional_mvout ! scalar
  integer(c_int) :: ts_do_additional_mvout_f_byte_size

  type(c_ptr) :: ts_do_mvout_only ! scalar
  integer(c_int) :: ts_do_mvout_only_f_byte_size

  type(c_ptr) :: ts_do_thermal_transonly ! scalar
  integer(c_int) :: ts_do_thermal_transonly_f_byte_size

  type(c_ptr) :: ts_do_observation_geometry ! scalar
  integer(c_int) :: ts_do_observation_geometry_f_byte_size

  type(c_ptr) :: ts_do_doublet_geometry ! scalar
  integer(c_int) :: ts_do_doublet_geometry_f_byte_size

  
end type lidort_modified_boolean_c

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_control_c
  type(c_ptr) :: ts_nmoments_input ! scalar
  integer(c_int) :: ts_nmoments_input_f_byte_size

  
end type lidort_modified_control_c

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_sunrays_c
  type(c_ptr) :: ts_nbeams ! scalar
  integer(c_int) :: ts_nbeams_f_byte_size

  type(c_ptr) :: ts_beam_szas ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_beam_szas_f_shapes
  integer(c_int) :: ts_beam_szas_f_byte_size

  
end type lidort_modified_sunrays_c

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_uservalues_c
  type(c_ptr) :: ts_n_user_relazms ! scalar
  integer(c_int) :: ts_n_user_relazms_f_byte_size

  type(c_ptr) :: ts_user_relazms ! dimension(MAX_USER_RELAZMS)
  integer(c_int), dimension(1) :: ts_user_relazms_f_shapes
  integer(c_int) :: ts_user_relazms_f_byte_size

  type(c_ptr) :: ts_n_user_streams ! scalar
  integer(c_int) :: ts_n_user_streams_f_byte_size

  type(c_ptr) :: ts_user_angles_input ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_user_angles_input_f_shapes
  integer(c_int) :: ts_user_angles_input_f_byte_size

  type(c_ptr) :: ts_user_levels ! dimension(MAX_USER_LEVELS)
  integer(c_int), dimension(1) :: ts_user_levels_f_shapes
  integer(c_int) :: ts_user_levels_f_byte_size

  type(c_ptr) :: ts_geometry_specheight ! scalar
  integer(c_int) :: ts_geometry_specheight_f_byte_size

  type(c_ptr) :: ts_n_user_obsgeoms ! scalar
  integer(c_int) :: ts_n_user_obsgeoms_f_byte_size

  type(c_ptr) :: ts_user_obsgeoms_input ! dimension(MAX_USER_OBSGEOMS, 3)
  integer(c_int), dimension(2) :: ts_user_obsgeoms_input_f_shapes
  integer(c_int) :: ts_user_obsgeoms_input_f_byte_size

  type(c_ptr) :: ts_n_user_doublets ! scalar
  integer(c_int) :: ts_n_user_doublets_f_byte_size

  type(c_ptr) :: ts_user_doublets ! dimension(MAX_USER_STREAMS, 2)
  integer(c_int), dimension(2) :: ts_user_doublets_f_shapes
  integer(c_int) :: ts_user_doublets_f_byte_size

  
end type lidort_modified_uservalues_c

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_chapman_c
  type(c_ptr) :: ts_earth_radius ! scalar
  integer(c_int) :: ts_earth_radius_f_byte_size

  
end type lidort_modified_chapman_c

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_optical_c
  type(c_ptr) :: ts_omega_total_input ! dimension(MAXLAYERS)
  integer(c_int), dimension(1) :: ts_omega_total_input_f_shapes
  integer(c_int) :: ts_omega_total_input_f_byte_size

  
end type lidort_modified_optical_c

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_inputs_c
  type(c_ptr) :: mbool ! scalar
  integer(c_int) :: mbool_f_byte_size

  type(c_ptr) :: mcont ! scalar
  integer(c_int) :: mcont_f_byte_size

  type(c_ptr) :: msunrays ! scalar
  integer(c_int) :: msunrays_f_byte_size

  type(c_ptr) :: muserval ! scalar
  integer(c_int) :: muserval_f_byte_size

  type(c_ptr) :: mchapman ! scalar
  integer(c_int) :: mchapman_f_byte_size

  type(c_ptr) :: moptical ! scalar
  integer(c_int) :: moptical_f_byte_size

  
end type lidort_modified_inputs_c


contains


subroutine set_lidort_pars(pars_struct) bind(C)
  type(lidort_pars_m_c) :: pars_struct
  integer :: len_idx

  do len_idx = 1, 5
    pars_struct%lidort_version_number(len_idx:len_idx) = LIDORT_VERSION_NUMBER(len_idx:len_idx)
  end do
  pars_struct%lidort_inunit = LIDORT_INUNIT
  pars_struct%lidort_scenunit = LIDORT_SCENUNIT
  pars_struct%lidort_funit = LIDORT_FUNIT
  pars_struct%lidort_resunit = LIDORT_RESUNIT
  pars_struct%lidort_errunit = LIDORT_ERRUNIT
  pars_struct%lidort_dbgunit = LIDORT_DBGUNIT
  pars_struct%max_messages = MAX_MESSAGES
  pars_struct%maxstreams = MAXSTREAMS
  pars_struct%maxlayers = MAXLAYERS
  pars_struct%maxfinelayers = MAXFINELAYERS
  pars_struct%maxmoments_input = MAXMOMENTS_INPUT
  pars_struct%max_thermal_coeffs = MAX_THERMAL_COEFFS
  pars_struct%maxbeams = MAXBEAMS
  pars_struct%max_user_streams = MAX_USER_STREAMS
  pars_struct%max_user_relazms = MAX_USER_RELAZMS
  pars_struct%max_user_obsgeoms = MAX_USER_OBSGEOMS
  pars_struct%max_user_levels = MAX_USER_LEVELS
  pars_struct%max_partlayers = MAX_PARTLAYERS
  pars_struct%max_taylor_terms = MAX_TAYLOR_TERMS
  pars_struct%max_directions = MAX_DIRECTIONS
  pars_struct%max_brdf_kernels = MAX_BRDF_KERNELS
  pars_struct%max_brdf_parameters = MAX_BRDF_PARAMETERS
  pars_struct%maxstreams_brdf = MAXSTREAMS_BRDF
  pars_struct%max_msrs_muquad = MAX_MSRS_MUQUAD
  pars_struct%max_msrs_phiquad = MAX_MSRS_PHIQUAD
  pars_struct%maxstreams_scaling = MAXSTREAMS_SCALING
  pars_struct%max_atmoswfs = MAX_ATMOSWFS
  pars_struct%max_surfacewfs = MAX_SURFACEWFS
  pars_struct%max_sleavewfs = MAX_SLEAVEWFS
  pars_struct%max_geometries = MAX_GEOMETRIES
  pars_struct%max_allstrms = MAX_ALLSTRMS
  pars_struct%max_allstrms_p1 = MAX_ALLSTRMS_P1
  pars_struct%maxmoments = MAXMOMENTS
  pars_struct%maxfourier = MAXFOURIER
  pars_struct%maxsthalf_brdf = MAXSTHALF_BRDF
  pars_struct%maxstreams_2 = MAXSTREAMS_2
  pars_struct%maxstreams_p1 = MAXSTREAMS_P1
  pars_struct%maxtotal = MAXTOTAL
  pars_struct%maxbandtotal = MAXBANDTOTAL
  pars_struct%one = ONE
  pars_struct%zero = ZERO
  pars_struct%onep5 = ONEP5
  pars_struct%two = TWO
  pars_struct%three = THREE
  pars_struct%four = FOUR
  pars_struct%quarter = QUARTER
  pars_struct%half = HALF
  pars_struct%minus_one = MINUS_ONE
  pars_struct%minus_two = MINUS_TWO
  pars_struct%pie = PIE
  pars_struct%deg_to_rad = DEG_TO_RAD
  pars_struct%pi2 = PI2
  pars_struct%pi4 = PI4
  pars_struct%pio2 = PIO2
  pars_struct%pio4 = PIO4
  pars_struct%eps3 = EPS3
  pars_struct%eps4 = EPS4
  pars_struct%eps5 = EPS5
  pars_struct%taylor_small = TAYLOR_SMALL
  pars_struct%smallnum = SMALLNUM
  pars_struct%bigexp = BIGEXP
  pars_struct%hopital_tolerance = HOPITAL_TOLERANCE
  pars_struct%omega_smallnum = OMEGA_SMALLNUM
  pars_struct%max_tau_spath = MAX_TAU_SPATH
  pars_struct%max_tau_upath = MAX_TAU_UPATH
  pars_struct%max_tau_qpath = MAX_TAU_QPATH
  pars_struct%lidort_serious = LIDORT_SERIOUS
  pars_struct%lidort_warning = LIDORT_WARNING
  pars_struct%lidort_info = LIDORT_INFO
  pars_struct%lidort_debug = LIDORT_DEBUG
  pars_struct%lidort_success = LIDORT_SUCCESS
  pars_struct%upidx = UPIDX
  pars_struct%dnidx = DNIDX
  pars_struct%lambertian_idx = LAMBERTIAN_IDX
  pars_struct%rossthin_idx = ROSSTHIN_IDX
  pars_struct%rossthick_idx = ROSSTHICK_IDX
  pars_struct%lisparse_idx = LISPARSE_IDX
  pars_struct%lidense_idx = LIDENSE_IDX
  pars_struct%hapke_idx = HAPKE_IDX
  pars_struct%roujean_idx = ROUJEAN_IDX
  pars_struct%rahman_idx = RAHMAN_IDX
  pars_struct%coxmunk_idx = COXMUNK_IDX
  pars_struct%bpdfsoil_idx = BPDFSOIL_IDX
  pars_struct%bpdfvegn_idx = BPDFVEGN_IDX
  pars_struct%bpdfndvi_idx = BPDFNDVI_IDX
  pars_struct%newcmglint_idx = NEWCMGLINT_IDX
  pars_struct%rtkhotspot_idx = RTKHOTSPOT_IDX
  pars_struct%modfresnel_idx = MODFRESNEL_IDX
  pars_struct%snowbrdf_idx = SNOWBRDF_IDX
  pars_struct%maxbrdf_idx = MAXBRDF_IDX
  
end subroutine set_lidort_pars



! Links to type: "brdf_linsup_inputs" from module: "brdf_lin_sup_inputs_def_m" in file: "brdf_lin_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_lin_sup_inputs_def_m, only : brdf_linsup_inputs

  type(brdf_linsup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_linsup_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_linsup_inputs_c_alloc_init

! Links to type: "brdf_linsup_inputs" from module: "brdf_lin_sup_inputs_def_m" in file: "brdf_lin_sup_inputs_def.F90"
! Initializes only with no allocation
subroutine brdf_linsup_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_lin_sup_inputs_def_m
  use lidort_pars_m

  type(brdf_linsup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  logical(kind=4), dimension(:), pointer :: bs_do_kernel_factor_wfs_ptr
  logical(kind=4), dimension(:,:), pointer :: bs_do_kernel_params_wfs_ptr
  logical(kind=4), dimension(:), pointer :: bs_do_kparams_derivs_ptr
  integer(c_int), pointer :: bs_n_surface_wfs_ptr
  integer(c_int), pointer :: bs_n_kernel_factor_wfs_ptr
  integer(c_int), pointer :: bs_n_kernel_params_wfs_ptr
  logical(kind=4), pointer :: bs_do_bsavalue_wf_ptr
  logical(kind=4), pointer :: bs_do_wsavalue_wf_ptr
  logical(kind=4), pointer :: bs_do_windspeed_wf_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_do_kernel_factor_wfs = .FALSE.
  bs_do_kernel_factor_wfs_ptr => fortran_type_f%bs_do_kernel_factor_wfs
  transfer_struct_c%bs_do_kernel_factor_wfs = c_loc(bs_do_kernel_factor_wfs_ptr(&
    lbound(fortran_type_f%bs_do_kernel_factor_wfs,1)))
  inquire(iolength=transfer_struct_c%bs_do_kernel_factor_wfs_f_byte_size) fortran_type_f%bs_do_kernel_factor_wfs(&
    lbound(fortran_type_f%bs_do_kernel_factor_wfs,1))
#ifdef ifort
  transfer_struct_c%bs_do_kernel_factor_wfs_f_byte_size = transfer_struct_c%bs_do_kernel_factor_wfs_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_do_kernel_factor_wfs_f_shapes(1) = size(fortran_type_f%bs_do_kernel_factor_wfs, 1)
  
  
  fortran_type_f%bs_do_kernel_params_wfs = .FALSE.
  bs_do_kernel_params_wfs_ptr => fortran_type_f%bs_do_kernel_params_wfs
  transfer_struct_c%bs_do_kernel_params_wfs = c_loc(bs_do_kernel_params_wfs_ptr(&
    lbound(fortran_type_f%bs_do_kernel_params_wfs,1),&
    lbound(fortran_type_f%bs_do_kernel_params_wfs,2)))
  inquire(iolength=transfer_struct_c%bs_do_kernel_params_wfs_f_byte_size) fortran_type_f%bs_do_kernel_params_wfs(&
    lbound(fortran_type_f%bs_do_kernel_params_wfs,1),&
    lbound(fortran_type_f%bs_do_kernel_params_wfs,2))
#ifdef ifort
  transfer_struct_c%bs_do_kernel_params_wfs_f_byte_size = transfer_struct_c%bs_do_kernel_params_wfs_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_do_kernel_params_wfs_f_shapes(1) = size(fortran_type_f%bs_do_kernel_params_wfs, 1)
  transfer_struct_c%bs_do_kernel_params_wfs_f_shapes(2) = size(fortran_type_f%bs_do_kernel_params_wfs, 2)
  
  
  fortran_type_f%bs_do_kparams_derivs = .FALSE.
  bs_do_kparams_derivs_ptr => fortran_type_f%bs_do_kparams_derivs
  transfer_struct_c%bs_do_kparams_derivs = c_loc(bs_do_kparams_derivs_ptr(&
    lbound(fortran_type_f%bs_do_kparams_derivs,1)))
  inquire(iolength=transfer_struct_c%bs_do_kparams_derivs_f_byte_size) fortran_type_f%bs_do_kparams_derivs(&
    lbound(fortran_type_f%bs_do_kparams_derivs,1))
#ifdef ifort
  transfer_struct_c%bs_do_kparams_derivs_f_byte_size = transfer_struct_c%bs_do_kparams_derivs_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_do_kparams_derivs_f_shapes(1) = size(fortran_type_f%bs_do_kparams_derivs, 1)
  
  
  fortran_type_f%bs_n_surface_wfs = 0
  bs_n_surface_wfs_ptr => fortran_type_f%bs_n_surface_wfs
  transfer_struct_c%bs_n_surface_wfs = c_loc(bs_n_surface_wfs_ptr)
  inquire(iolength=transfer_struct_c%bs_n_surface_wfs_f_byte_size) fortran_type_f%bs_n_surface_wfs
#ifdef ifort
  transfer_struct_c%bs_n_surface_wfs_f_byte_size = transfer_struct_c%bs_n_surface_wfs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_n_kernel_factor_wfs = 0
  bs_n_kernel_factor_wfs_ptr => fortran_type_f%bs_n_kernel_factor_wfs
  transfer_struct_c%bs_n_kernel_factor_wfs = c_loc(bs_n_kernel_factor_wfs_ptr)
  inquire(iolength=transfer_struct_c%bs_n_kernel_factor_wfs_f_byte_size) fortran_type_f%bs_n_kernel_factor_wfs
#ifdef ifort
  transfer_struct_c%bs_n_kernel_factor_wfs_f_byte_size = transfer_struct_c%bs_n_kernel_factor_wfs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_n_kernel_params_wfs = 0
  bs_n_kernel_params_wfs_ptr => fortran_type_f%bs_n_kernel_params_wfs
  transfer_struct_c%bs_n_kernel_params_wfs = c_loc(bs_n_kernel_params_wfs_ptr)
  inquire(iolength=transfer_struct_c%bs_n_kernel_params_wfs_f_byte_size) fortran_type_f%bs_n_kernel_params_wfs
#ifdef ifort
  transfer_struct_c%bs_n_kernel_params_wfs_f_byte_size = transfer_struct_c%bs_n_kernel_params_wfs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_bsavalue_wf = .FALSE.
  bs_do_bsavalue_wf_ptr => fortran_type_f%bs_do_bsavalue_wf
  transfer_struct_c%bs_do_bsavalue_wf = c_loc(bs_do_bsavalue_wf_ptr)
  inquire(iolength=transfer_struct_c%bs_do_bsavalue_wf_f_byte_size) fortran_type_f%bs_do_bsavalue_wf
#ifdef ifort
  transfer_struct_c%bs_do_bsavalue_wf_f_byte_size = transfer_struct_c%bs_do_bsavalue_wf_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_wsavalue_wf = .FALSE.
  bs_do_wsavalue_wf_ptr => fortran_type_f%bs_do_wsavalue_wf
  transfer_struct_c%bs_do_wsavalue_wf = c_loc(bs_do_wsavalue_wf_ptr)
  inquire(iolength=transfer_struct_c%bs_do_wsavalue_wf_f_byte_size) fortran_type_f%bs_do_wsavalue_wf
#ifdef ifort
  transfer_struct_c%bs_do_wsavalue_wf_f_byte_size = transfer_struct_c%bs_do_wsavalue_wf_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_windspeed_wf = .FALSE.
  bs_do_windspeed_wf_ptr => fortran_type_f%bs_do_windspeed_wf
  transfer_struct_c%bs_do_windspeed_wf = c_loc(bs_do_windspeed_wf_ptr)
  inquire(iolength=transfer_struct_c%bs_do_windspeed_wf_f_byte_size) fortran_type_f%bs_do_windspeed_wf
#ifdef ifort
  transfer_struct_c%bs_do_windspeed_wf_f_byte_size = transfer_struct_c%bs_do_windspeed_wf_f_byte_size * 4
#endif
  
  
  
end subroutine brdf_linsup_inputs_c_init_only

subroutine brdf_linsup_inputs_c_destroy(fortran_type_c) bind(C)
  use brdf_lin_sup_inputs_def_m, only : brdf_linsup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_linsup_inputs_c_destroy

subroutine brdf_linsup_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_lin_sup_inputs_def_m, only : brdf_linsup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_linsup_inputs), pointer :: fortran_type_f_from
  type(brdf_linsup_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_do_kernel_factor_wfs = fortran_type_f_from%bs_do_kernel_factor_wfs
  fortran_type_f_to%bs_do_kernel_params_wfs = fortran_type_f_from%bs_do_kernel_params_wfs
  fortran_type_f_to%bs_do_kparams_derivs = fortran_type_f_from%bs_do_kparams_derivs
  fortran_type_f_to%bs_n_surface_wfs = fortran_type_f_from%bs_n_surface_wfs
  fortran_type_f_to%bs_n_kernel_factor_wfs = fortran_type_f_from%bs_n_kernel_factor_wfs
  fortran_type_f_to%bs_n_kernel_params_wfs = fortran_type_f_from%bs_n_kernel_params_wfs
  fortran_type_f_to%bs_do_bsavalue_wf = fortran_type_f_from%bs_do_bsavalue_wf
  fortran_type_f_to%bs_do_wsavalue_wf = fortran_type_f_from%bs_do_wsavalue_wf
  fortran_type_f_to%bs_do_windspeed_wf = fortran_type_f_from%bs_do_windspeed_wf
  

end subroutine brdf_linsup_inputs_c_copy

! Links to type: "brdf_linsup_outputs" from module: "brdf_lin_sup_outputs_def_m" in file: "brdf_lin_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_lin_sup_outputs_def_m, only : brdf_linsup_outputs

  type(brdf_linsup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_linsup_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_linsup_outputs_c_alloc_init

! Links to type: "brdf_linsup_outputs" from module: "brdf_lin_sup_outputs_def_m" in file: "brdf_lin_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_linsup_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_lin_sup_outputs_def_m
  use lidort_pars_m

  type(brdf_linsup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_dbounce_brdfunc_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_brdf_f_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_user_brdf_f_ptr
  real(c_double), dimension(:,:), pointer :: bs_ls_emissivity_ptr
  real(c_double), dimension(:,:), pointer :: bs_ls_user_emissivity_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_ls_dbounce_brdfunc = 0_fpk
  bs_ls_dbounce_brdfunc_ptr => fortran_type_f%bs_ls_dbounce_brdfunc
  transfer_struct_c%bs_ls_dbounce_brdfunc = c_loc(bs_ls_dbounce_brdfunc_ptr(&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,1),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,2),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,3),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,4)))
  inquire(iolength=transfer_struct_c%bs_ls_dbounce_brdfunc_f_byte_size) fortran_type_f%bs_ls_dbounce_brdfunc(&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,1),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,2),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,3),&
    lbound(fortran_type_f%bs_ls_dbounce_brdfunc,4))
#ifdef ifort
  transfer_struct_c%bs_ls_dbounce_brdfunc_f_byte_size = transfer_struct_c%bs_ls_dbounce_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_dbounce_brdfunc_f_shapes(1) = size(fortran_type_f%bs_ls_dbounce_brdfunc, 1)
  transfer_struct_c%bs_ls_dbounce_brdfunc_f_shapes(2) = size(fortran_type_f%bs_ls_dbounce_brdfunc, 2)
  transfer_struct_c%bs_ls_dbounce_brdfunc_f_shapes(3) = size(fortran_type_f%bs_ls_dbounce_brdfunc, 3)
  transfer_struct_c%bs_ls_dbounce_brdfunc_f_shapes(4) = size(fortran_type_f%bs_ls_dbounce_brdfunc, 4)
  
  
  fortran_type_f%bs_ls_brdf_f_0 = 0_fpk
  bs_ls_brdf_f_0_ptr => fortran_type_f%bs_ls_brdf_f_0
  transfer_struct_c%bs_ls_brdf_f_0 = c_loc(bs_ls_brdf_f_0_ptr(&
    lbound(fortran_type_f%bs_ls_brdf_f_0,1),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,2),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,3),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,4)))
  inquire(iolength=transfer_struct_c%bs_ls_brdf_f_0_f_byte_size) fortran_type_f%bs_ls_brdf_f_0(&
    lbound(fortran_type_f%bs_ls_brdf_f_0,1),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,2),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,3),&
    lbound(fortran_type_f%bs_ls_brdf_f_0,4))
#ifdef ifort
  transfer_struct_c%bs_ls_brdf_f_0_f_byte_size = transfer_struct_c%bs_ls_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_brdf_f_0_f_shapes(1) = size(fortran_type_f%bs_ls_brdf_f_0, 1)
  transfer_struct_c%bs_ls_brdf_f_0_f_shapes(2) = size(fortran_type_f%bs_ls_brdf_f_0, 2)
  transfer_struct_c%bs_ls_brdf_f_0_f_shapes(3) = size(fortran_type_f%bs_ls_brdf_f_0, 3)
  transfer_struct_c%bs_ls_brdf_f_0_f_shapes(4) = size(fortran_type_f%bs_ls_brdf_f_0, 4)
  
  
  fortran_type_f%bs_ls_brdf_f = 0_fpk
  bs_ls_brdf_f_ptr => fortran_type_f%bs_ls_brdf_f
  transfer_struct_c%bs_ls_brdf_f = c_loc(bs_ls_brdf_f_ptr(&
    lbound(fortran_type_f%bs_ls_brdf_f,1),&
    lbound(fortran_type_f%bs_ls_brdf_f,2),&
    lbound(fortran_type_f%bs_ls_brdf_f,3),&
    lbound(fortran_type_f%bs_ls_brdf_f,4)))
  inquire(iolength=transfer_struct_c%bs_ls_brdf_f_f_byte_size) fortran_type_f%bs_ls_brdf_f(&
    lbound(fortran_type_f%bs_ls_brdf_f,1),&
    lbound(fortran_type_f%bs_ls_brdf_f,2),&
    lbound(fortran_type_f%bs_ls_brdf_f,3),&
    lbound(fortran_type_f%bs_ls_brdf_f,4))
#ifdef ifort
  transfer_struct_c%bs_ls_brdf_f_f_byte_size = transfer_struct_c%bs_ls_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_brdf_f_f_shapes(1) = size(fortran_type_f%bs_ls_brdf_f, 1)
  transfer_struct_c%bs_ls_brdf_f_f_shapes(2) = size(fortran_type_f%bs_ls_brdf_f, 2)
  transfer_struct_c%bs_ls_brdf_f_f_shapes(3) = size(fortran_type_f%bs_ls_brdf_f, 3)
  transfer_struct_c%bs_ls_brdf_f_f_shapes(4) = size(fortran_type_f%bs_ls_brdf_f, 4)
  
  
  fortran_type_f%bs_ls_user_brdf_f_0 = 0_fpk
  bs_ls_user_brdf_f_0_ptr => fortran_type_f%bs_ls_user_brdf_f_0
  transfer_struct_c%bs_ls_user_brdf_f_0 = c_loc(bs_ls_user_brdf_f_0_ptr(&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,1),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,2),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,3),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,4)))
  inquire(iolength=transfer_struct_c%bs_ls_user_brdf_f_0_f_byte_size) fortran_type_f%bs_ls_user_brdf_f_0(&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,1),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,2),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,3),&
    lbound(fortran_type_f%bs_ls_user_brdf_f_0,4))
#ifdef ifort
  transfer_struct_c%bs_ls_user_brdf_f_0_f_byte_size = transfer_struct_c%bs_ls_user_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_user_brdf_f_0_f_shapes(1) = size(fortran_type_f%bs_ls_user_brdf_f_0, 1)
  transfer_struct_c%bs_ls_user_brdf_f_0_f_shapes(2) = size(fortran_type_f%bs_ls_user_brdf_f_0, 2)
  transfer_struct_c%bs_ls_user_brdf_f_0_f_shapes(3) = size(fortran_type_f%bs_ls_user_brdf_f_0, 3)
  transfer_struct_c%bs_ls_user_brdf_f_0_f_shapes(4) = size(fortran_type_f%bs_ls_user_brdf_f_0, 4)
  
  
  fortran_type_f%bs_ls_user_brdf_f = 0_fpk
  bs_ls_user_brdf_f_ptr => fortran_type_f%bs_ls_user_brdf_f
  transfer_struct_c%bs_ls_user_brdf_f = c_loc(bs_ls_user_brdf_f_ptr(&
    lbound(fortran_type_f%bs_ls_user_brdf_f,1),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,2),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,3),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,4)))
  inquire(iolength=transfer_struct_c%bs_ls_user_brdf_f_f_byte_size) fortran_type_f%bs_ls_user_brdf_f(&
    lbound(fortran_type_f%bs_ls_user_brdf_f,1),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,2),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,3),&
    lbound(fortran_type_f%bs_ls_user_brdf_f,4))
#ifdef ifort
  transfer_struct_c%bs_ls_user_brdf_f_f_byte_size = transfer_struct_c%bs_ls_user_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_user_brdf_f_f_shapes(1) = size(fortran_type_f%bs_ls_user_brdf_f, 1)
  transfer_struct_c%bs_ls_user_brdf_f_f_shapes(2) = size(fortran_type_f%bs_ls_user_brdf_f, 2)
  transfer_struct_c%bs_ls_user_brdf_f_f_shapes(3) = size(fortran_type_f%bs_ls_user_brdf_f, 3)
  transfer_struct_c%bs_ls_user_brdf_f_f_shapes(4) = size(fortran_type_f%bs_ls_user_brdf_f, 4)
  
  
  fortran_type_f%bs_ls_emissivity = 0_fpk
  bs_ls_emissivity_ptr => fortran_type_f%bs_ls_emissivity
  transfer_struct_c%bs_ls_emissivity = c_loc(bs_ls_emissivity_ptr(&
    lbound(fortran_type_f%bs_ls_emissivity,1),&
    lbound(fortran_type_f%bs_ls_emissivity,2)))
  inquire(iolength=transfer_struct_c%bs_ls_emissivity_f_byte_size) fortran_type_f%bs_ls_emissivity(&
    lbound(fortran_type_f%bs_ls_emissivity,1),&
    lbound(fortran_type_f%bs_ls_emissivity,2))
#ifdef ifort
  transfer_struct_c%bs_ls_emissivity_f_byte_size = transfer_struct_c%bs_ls_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_emissivity_f_shapes(1) = size(fortran_type_f%bs_ls_emissivity, 1)
  transfer_struct_c%bs_ls_emissivity_f_shapes(2) = size(fortran_type_f%bs_ls_emissivity, 2)
  
  
  fortran_type_f%bs_ls_user_emissivity = 0_fpk
  bs_ls_user_emissivity_ptr => fortran_type_f%bs_ls_user_emissivity
  transfer_struct_c%bs_ls_user_emissivity = c_loc(bs_ls_user_emissivity_ptr(&
    lbound(fortran_type_f%bs_ls_user_emissivity,1),&
    lbound(fortran_type_f%bs_ls_user_emissivity,2)))
  inquire(iolength=transfer_struct_c%bs_ls_user_emissivity_f_byte_size) fortran_type_f%bs_ls_user_emissivity(&
    lbound(fortran_type_f%bs_ls_user_emissivity,1),&
    lbound(fortran_type_f%bs_ls_user_emissivity,2))
#ifdef ifort
  transfer_struct_c%bs_ls_user_emissivity_f_byte_size = transfer_struct_c%bs_ls_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_user_emissivity_f_shapes(1) = size(fortran_type_f%bs_ls_user_emissivity, 1)
  transfer_struct_c%bs_ls_user_emissivity_f_shapes(2) = size(fortran_type_f%bs_ls_user_emissivity, 2)
  
  
end subroutine brdf_linsup_outputs_c_init_only

subroutine brdf_linsup_outputs_c_destroy(fortran_type_c) bind(C)
  use brdf_lin_sup_outputs_def_m, only : brdf_linsup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_linsup_outputs_c_destroy

subroutine brdf_linsup_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_lin_sup_outputs_def_m, only : brdf_linsup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_linsup_outputs), pointer :: fortran_type_f_from
  type(brdf_linsup_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_ls_dbounce_brdfunc = fortran_type_f_from%bs_ls_dbounce_brdfunc
  fortran_type_f_to%bs_ls_brdf_f_0 = fortran_type_f_from%bs_ls_brdf_f_0
  fortran_type_f_to%bs_ls_brdf_f = fortran_type_f_from%bs_ls_brdf_f
  fortran_type_f_to%bs_ls_user_brdf_f_0 = fortran_type_f_from%bs_ls_user_brdf_f_0
  fortran_type_f_to%bs_ls_user_brdf_f = fortran_type_f_from%bs_ls_user_brdf_f
  fortran_type_f_to%bs_ls_emissivity = fortran_type_f_from%bs_ls_emissivity
  fortran_type_f_to%bs_ls_user_emissivity = fortran_type_f_from%bs_ls_user_emissivity
  

end subroutine brdf_linsup_outputs_c_copy

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_inputs_def_m, only : brdf_sup_inputs

  type(brdf_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_sup_inputs_c_alloc_init

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.F90"
! Initializes only with no allocation
subroutine brdf_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_inputs_def_m
  use lidort_pars_m

  type(brdf_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  logical(kind=4), pointer :: bs_do_brdf_surface_ptr
  logical(kind=4), pointer :: bs_do_surface_emission_ptr
  logical(kind=4), pointer :: bs_do_solar_sources_ptr
  logical(kind=4), pointer :: bs_do_user_streams_ptr
  logical(kind=4), pointer :: bs_do_user_obsgeoms_ptr
  logical(kind=4), pointer :: bs_do_doublet_geometry_ptr
  integer(c_int), pointer :: bs_nstreams_ptr
  integer(c_int), pointer :: bs_nbeams_ptr
  real(c_double), dimension(:), pointer :: bs_beam_szas_ptr
  integer(c_int), pointer :: bs_n_user_relazms_ptr
  real(c_double), dimension(:), pointer :: bs_user_relazms_ptr
  integer(c_int), pointer :: bs_n_user_streams_ptr
  real(c_double), dimension(:), pointer :: bs_user_angles_input_ptr
  integer(c_int), pointer :: bs_n_user_obsgeoms_ptr
  real(c_double), dimension(:,:), pointer :: bs_user_obsgeoms_ptr
  integer(c_int), pointer :: bs_n_user_doublets_ptr
  real(c_double), dimension(:,:), pointer :: bs_user_doublets_ptr
  integer(c_int), pointer :: bs_n_brdf_kernels_ptr
  integer(c_int), dimension(:), pointer :: bs_which_brdf_ptr
  integer(c_int), dimension(:), pointer :: bs_n_brdf_parameters_ptr
  real(c_double), dimension(:,:), pointer :: bs_brdf_parameters_ptr
  logical(kind=4), dimension(:), pointer :: bs_lambertian_kernel_flag_ptr
  real(c_double), dimension(:), pointer :: bs_brdf_factors_ptr
  integer(c_int), pointer :: bs_nstreams_brdf_ptr
  logical(kind=4), pointer :: bs_do_shadow_effect_ptr
  logical(kind=4), pointer :: bs_do_directbounce_only_ptr
  logical(kind=4), pointer :: bs_do_wsabsa_output_ptr
  logical(kind=4), pointer :: bs_do_wsa_scaling_ptr
  logical(kind=4), pointer :: bs_do_bsa_scaling_ptr
  real(c_double), pointer :: bs_wsa_value_ptr
  real(c_double), pointer :: bs_bsa_value_ptr
  logical(kind=4), pointer :: bs_do_newcmglint_ptr
  real(c_double), pointer :: bs_salinity_ptr
  real(c_double), pointer :: bs_wavelength_ptr
  real(c_double), pointer :: bs_windspeed_ptr
  real(c_double), dimension(:), pointer :: bs_winddir_ptr
  logical(kind=4), pointer :: bs_do_glintshadow_ptr
  logical(kind=4), pointer :: bs_do_foamoption_ptr
  logical(kind=4), pointer :: bs_do_facetisotropy_ptr
  logical(kind=4), pointer :: bs_do_glitter_msrcorr_ptr
  logical(kind=4), pointer :: bs_do_glitter_msrcorr_dbonly_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_order_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_nmuquad_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_nphiquad_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_do_brdf_surface = .FALSE.
  bs_do_brdf_surface_ptr => fortran_type_f%bs_do_brdf_surface
  transfer_struct_c%bs_do_brdf_surface = c_loc(bs_do_brdf_surface_ptr)
  inquire(iolength=transfer_struct_c%bs_do_brdf_surface_f_byte_size) fortran_type_f%bs_do_brdf_surface
#ifdef ifort
  transfer_struct_c%bs_do_brdf_surface_f_byte_size = transfer_struct_c%bs_do_brdf_surface_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_surface_emission = .FALSE.
  bs_do_surface_emission_ptr => fortran_type_f%bs_do_surface_emission
  transfer_struct_c%bs_do_surface_emission = c_loc(bs_do_surface_emission_ptr)
  inquire(iolength=transfer_struct_c%bs_do_surface_emission_f_byte_size) fortran_type_f%bs_do_surface_emission
#ifdef ifort
  transfer_struct_c%bs_do_surface_emission_f_byte_size = transfer_struct_c%bs_do_surface_emission_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_solar_sources = .FALSE.
  bs_do_solar_sources_ptr => fortran_type_f%bs_do_solar_sources
  transfer_struct_c%bs_do_solar_sources = c_loc(bs_do_solar_sources_ptr)
  inquire(iolength=transfer_struct_c%bs_do_solar_sources_f_byte_size) fortran_type_f%bs_do_solar_sources
#ifdef ifort
  transfer_struct_c%bs_do_solar_sources_f_byte_size = transfer_struct_c%bs_do_solar_sources_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_user_streams = .FALSE.
  bs_do_user_streams_ptr => fortran_type_f%bs_do_user_streams
  transfer_struct_c%bs_do_user_streams = c_loc(bs_do_user_streams_ptr)
  inquire(iolength=transfer_struct_c%bs_do_user_streams_f_byte_size) fortran_type_f%bs_do_user_streams
#ifdef ifort
  transfer_struct_c%bs_do_user_streams_f_byte_size = transfer_struct_c%bs_do_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_user_obsgeoms = .FALSE.
  bs_do_user_obsgeoms_ptr => fortran_type_f%bs_do_user_obsgeoms
  transfer_struct_c%bs_do_user_obsgeoms = c_loc(bs_do_user_obsgeoms_ptr)
  inquire(iolength=transfer_struct_c%bs_do_user_obsgeoms_f_byte_size) fortran_type_f%bs_do_user_obsgeoms
#ifdef ifort
  transfer_struct_c%bs_do_user_obsgeoms_f_byte_size = transfer_struct_c%bs_do_user_obsgeoms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_doublet_geometry = .FALSE.
  bs_do_doublet_geometry_ptr => fortran_type_f%bs_do_doublet_geometry
  transfer_struct_c%bs_do_doublet_geometry = c_loc(bs_do_doublet_geometry_ptr)
  inquire(iolength=transfer_struct_c%bs_do_doublet_geometry_f_byte_size) fortran_type_f%bs_do_doublet_geometry
#ifdef ifort
  transfer_struct_c%bs_do_doublet_geometry_f_byte_size = transfer_struct_c%bs_do_doublet_geometry_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_nstreams = 0
  bs_nstreams_ptr => fortran_type_f%bs_nstreams
  transfer_struct_c%bs_nstreams = c_loc(bs_nstreams_ptr)
  inquire(iolength=transfer_struct_c%bs_nstreams_f_byte_size) fortran_type_f%bs_nstreams
#ifdef ifort
  transfer_struct_c%bs_nstreams_f_byte_size = transfer_struct_c%bs_nstreams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_nbeams = 0
  bs_nbeams_ptr => fortran_type_f%bs_nbeams
  transfer_struct_c%bs_nbeams = c_loc(bs_nbeams_ptr)
  inquire(iolength=transfer_struct_c%bs_nbeams_f_byte_size) fortran_type_f%bs_nbeams
#ifdef ifort
  transfer_struct_c%bs_nbeams_f_byte_size = transfer_struct_c%bs_nbeams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_beam_szas = 0_fpk
  bs_beam_szas_ptr => fortran_type_f%bs_beam_szas
  transfer_struct_c%bs_beam_szas = c_loc(bs_beam_szas_ptr(&
    lbound(fortran_type_f%bs_beam_szas,1)))
  inquire(iolength=transfer_struct_c%bs_beam_szas_f_byte_size) fortran_type_f%bs_beam_szas(&
    lbound(fortran_type_f%bs_beam_szas,1))
#ifdef ifort
  transfer_struct_c%bs_beam_szas_f_byte_size = transfer_struct_c%bs_beam_szas_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_beam_szas_f_shapes(1) = size(fortran_type_f%bs_beam_szas, 1)
  
  
  fortran_type_f%bs_n_user_relazms = 0
  bs_n_user_relazms_ptr => fortran_type_f%bs_n_user_relazms
  transfer_struct_c%bs_n_user_relazms = c_loc(bs_n_user_relazms_ptr)
  inquire(iolength=transfer_struct_c%bs_n_user_relazms_f_byte_size) fortran_type_f%bs_n_user_relazms
#ifdef ifort
  transfer_struct_c%bs_n_user_relazms_f_byte_size = transfer_struct_c%bs_n_user_relazms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_user_relazms = 0_fpk
  bs_user_relazms_ptr => fortran_type_f%bs_user_relazms
  transfer_struct_c%bs_user_relazms = c_loc(bs_user_relazms_ptr(&
    lbound(fortran_type_f%bs_user_relazms,1)))
  inquire(iolength=transfer_struct_c%bs_user_relazms_f_byte_size) fortran_type_f%bs_user_relazms(&
    lbound(fortran_type_f%bs_user_relazms,1))
#ifdef ifort
  transfer_struct_c%bs_user_relazms_f_byte_size = transfer_struct_c%bs_user_relazms_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_relazms_f_shapes(1) = size(fortran_type_f%bs_user_relazms, 1)
  
  
  fortran_type_f%bs_n_user_streams = 0
  bs_n_user_streams_ptr => fortran_type_f%bs_n_user_streams
  transfer_struct_c%bs_n_user_streams = c_loc(bs_n_user_streams_ptr)
  inquire(iolength=transfer_struct_c%bs_n_user_streams_f_byte_size) fortran_type_f%bs_n_user_streams
#ifdef ifort
  transfer_struct_c%bs_n_user_streams_f_byte_size = transfer_struct_c%bs_n_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_user_angles_input = 0_fpk
  bs_user_angles_input_ptr => fortran_type_f%bs_user_angles_input
  transfer_struct_c%bs_user_angles_input = c_loc(bs_user_angles_input_ptr(&
    lbound(fortran_type_f%bs_user_angles_input,1)))
  inquire(iolength=transfer_struct_c%bs_user_angles_input_f_byte_size) fortran_type_f%bs_user_angles_input(&
    lbound(fortran_type_f%bs_user_angles_input,1))
#ifdef ifort
  transfer_struct_c%bs_user_angles_input_f_byte_size = transfer_struct_c%bs_user_angles_input_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_angles_input_f_shapes(1) = size(fortran_type_f%bs_user_angles_input, 1)
  
  
  fortran_type_f%bs_n_user_obsgeoms = 0
  bs_n_user_obsgeoms_ptr => fortran_type_f%bs_n_user_obsgeoms
  transfer_struct_c%bs_n_user_obsgeoms = c_loc(bs_n_user_obsgeoms_ptr)
  inquire(iolength=transfer_struct_c%bs_n_user_obsgeoms_f_byte_size) fortran_type_f%bs_n_user_obsgeoms
#ifdef ifort
  transfer_struct_c%bs_n_user_obsgeoms_f_byte_size = transfer_struct_c%bs_n_user_obsgeoms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_user_obsgeoms = 0_fpk
  bs_user_obsgeoms_ptr => fortran_type_f%bs_user_obsgeoms
  transfer_struct_c%bs_user_obsgeoms = c_loc(bs_user_obsgeoms_ptr(&
    lbound(fortran_type_f%bs_user_obsgeoms,1),&
    lbound(fortran_type_f%bs_user_obsgeoms,2)))
  inquire(iolength=transfer_struct_c%bs_user_obsgeoms_f_byte_size) fortran_type_f%bs_user_obsgeoms(&
    lbound(fortran_type_f%bs_user_obsgeoms,1),&
    lbound(fortran_type_f%bs_user_obsgeoms,2))
#ifdef ifort
  transfer_struct_c%bs_user_obsgeoms_f_byte_size = transfer_struct_c%bs_user_obsgeoms_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_obsgeoms_f_shapes(1) = size(fortran_type_f%bs_user_obsgeoms, 1)
  transfer_struct_c%bs_user_obsgeoms_f_shapes(2) = size(fortran_type_f%bs_user_obsgeoms, 2)
  
  
  fortran_type_f%bs_n_user_doublets = 0
  bs_n_user_doublets_ptr => fortran_type_f%bs_n_user_doublets
  transfer_struct_c%bs_n_user_doublets = c_loc(bs_n_user_doublets_ptr)
  inquire(iolength=transfer_struct_c%bs_n_user_doublets_f_byte_size) fortran_type_f%bs_n_user_doublets
#ifdef ifort
  transfer_struct_c%bs_n_user_doublets_f_byte_size = transfer_struct_c%bs_n_user_doublets_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_user_doublets = 0_fpk
  bs_user_doublets_ptr => fortran_type_f%bs_user_doublets
  transfer_struct_c%bs_user_doublets = c_loc(bs_user_doublets_ptr(&
    lbound(fortran_type_f%bs_user_doublets,1),&
    lbound(fortran_type_f%bs_user_doublets,2)))
  inquire(iolength=transfer_struct_c%bs_user_doublets_f_byte_size) fortran_type_f%bs_user_doublets(&
    lbound(fortran_type_f%bs_user_doublets,1),&
    lbound(fortran_type_f%bs_user_doublets,2))
#ifdef ifort
  transfer_struct_c%bs_user_doublets_f_byte_size = transfer_struct_c%bs_user_doublets_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_doublets_f_shapes(1) = size(fortran_type_f%bs_user_doublets, 1)
  transfer_struct_c%bs_user_doublets_f_shapes(2) = size(fortran_type_f%bs_user_doublets, 2)
  
  
  fortran_type_f%bs_n_brdf_kernels = 0
  bs_n_brdf_kernels_ptr => fortran_type_f%bs_n_brdf_kernels
  transfer_struct_c%bs_n_brdf_kernels = c_loc(bs_n_brdf_kernels_ptr)
  inquire(iolength=transfer_struct_c%bs_n_brdf_kernels_f_byte_size) fortran_type_f%bs_n_brdf_kernels
#ifdef ifort
  transfer_struct_c%bs_n_brdf_kernels_f_byte_size = transfer_struct_c%bs_n_brdf_kernels_f_byte_size * 4
#endif
  
  
  fortran_type_f%bs_brdf_names = ''
  transfer_struct_c%bs_brdf_names_f_len = len(fortran_type_f%bs_brdf_names)
  transfer_struct_c%bs_brdf_names_f_shapes(1) = size(fortran_type_f%bs_brdf_names, 1)
  
  
  fortran_type_f%bs_which_brdf = 0
  bs_which_brdf_ptr => fortran_type_f%bs_which_brdf
  transfer_struct_c%bs_which_brdf = c_loc(bs_which_brdf_ptr(&
    lbound(fortran_type_f%bs_which_brdf,1)))
  inquire(iolength=transfer_struct_c%bs_which_brdf_f_byte_size) fortran_type_f%bs_which_brdf(&
    lbound(fortran_type_f%bs_which_brdf,1))
#ifdef ifort
  transfer_struct_c%bs_which_brdf_f_byte_size = transfer_struct_c%bs_which_brdf_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_which_brdf_f_shapes(1) = size(fortran_type_f%bs_which_brdf, 1)
  
  
  fortran_type_f%bs_n_brdf_parameters = 0
  bs_n_brdf_parameters_ptr => fortran_type_f%bs_n_brdf_parameters
  transfer_struct_c%bs_n_brdf_parameters = c_loc(bs_n_brdf_parameters_ptr(&
    lbound(fortran_type_f%bs_n_brdf_parameters,1)))
  inquire(iolength=transfer_struct_c%bs_n_brdf_parameters_f_byte_size) fortran_type_f%bs_n_brdf_parameters(&
    lbound(fortran_type_f%bs_n_brdf_parameters,1))
#ifdef ifort
  transfer_struct_c%bs_n_brdf_parameters_f_byte_size = transfer_struct_c%bs_n_brdf_parameters_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_n_brdf_parameters_f_shapes(1) = size(fortran_type_f%bs_n_brdf_parameters, 1)
  
  
  fortran_type_f%bs_brdf_parameters = 0_fpk
  bs_brdf_parameters_ptr => fortran_type_f%bs_brdf_parameters
  transfer_struct_c%bs_brdf_parameters = c_loc(bs_brdf_parameters_ptr(&
    lbound(fortran_type_f%bs_brdf_parameters,1),&
    lbound(fortran_type_f%bs_brdf_parameters,2)))
  inquire(iolength=transfer_struct_c%bs_brdf_parameters_f_byte_size) fortran_type_f%bs_brdf_parameters(&
    lbound(fortran_type_f%bs_brdf_parameters,1),&
    lbound(fortran_type_f%bs_brdf_parameters,2))
#ifdef ifort
  transfer_struct_c%bs_brdf_parameters_f_byte_size = transfer_struct_c%bs_brdf_parameters_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_brdf_parameters_f_shapes(1) = size(fortran_type_f%bs_brdf_parameters, 1)
  transfer_struct_c%bs_brdf_parameters_f_shapes(2) = size(fortran_type_f%bs_brdf_parameters, 2)
  
  
  fortran_type_f%bs_lambertian_kernel_flag = .FALSE.
  bs_lambertian_kernel_flag_ptr => fortran_type_f%bs_lambertian_kernel_flag
  transfer_struct_c%bs_lambertian_kernel_flag = c_loc(bs_lambertian_kernel_flag_ptr(&
    lbound(fortran_type_f%bs_lambertian_kernel_flag,1)))
  inquire(iolength=transfer_struct_c%bs_lambertian_kernel_flag_f_byte_size) fortran_type_f%bs_lambertian_kernel_flag(&
    lbound(fortran_type_f%bs_lambertian_kernel_flag,1))
#ifdef ifort
  transfer_struct_c%bs_lambertian_kernel_flag_f_byte_size = transfer_struct_c%bs_lambertian_kernel_flag_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_lambertian_kernel_flag_f_shapes(1) = size(fortran_type_f%bs_lambertian_kernel_flag, 1)
  
  
  fortran_type_f%bs_brdf_factors = 0_fpk
  bs_brdf_factors_ptr => fortran_type_f%bs_brdf_factors
  transfer_struct_c%bs_brdf_factors = c_loc(bs_brdf_factors_ptr(&
    lbound(fortran_type_f%bs_brdf_factors,1)))
  inquire(iolength=transfer_struct_c%bs_brdf_factors_f_byte_size) fortran_type_f%bs_brdf_factors(&
    lbound(fortran_type_f%bs_brdf_factors,1))
#ifdef ifort
  transfer_struct_c%bs_brdf_factors_f_byte_size = transfer_struct_c%bs_brdf_factors_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_brdf_factors_f_shapes(1) = size(fortran_type_f%bs_brdf_factors, 1)
  
  
  fortran_type_f%bs_nstreams_brdf = 0
  bs_nstreams_brdf_ptr => fortran_type_f%bs_nstreams_brdf
  transfer_struct_c%bs_nstreams_brdf = c_loc(bs_nstreams_brdf_ptr)
  inquire(iolength=transfer_struct_c%bs_nstreams_brdf_f_byte_size) fortran_type_f%bs_nstreams_brdf
#ifdef ifort
  transfer_struct_c%bs_nstreams_brdf_f_byte_size = transfer_struct_c%bs_nstreams_brdf_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_shadow_effect = .FALSE.
  bs_do_shadow_effect_ptr => fortran_type_f%bs_do_shadow_effect
  transfer_struct_c%bs_do_shadow_effect = c_loc(bs_do_shadow_effect_ptr)
  inquire(iolength=transfer_struct_c%bs_do_shadow_effect_f_byte_size) fortran_type_f%bs_do_shadow_effect
#ifdef ifort
  transfer_struct_c%bs_do_shadow_effect_f_byte_size = transfer_struct_c%bs_do_shadow_effect_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_directbounce_only = .FALSE.
  bs_do_directbounce_only_ptr => fortran_type_f%bs_do_directbounce_only
  transfer_struct_c%bs_do_directbounce_only = c_loc(bs_do_directbounce_only_ptr)
  inquire(iolength=transfer_struct_c%bs_do_directbounce_only_f_byte_size) fortran_type_f%bs_do_directbounce_only
#ifdef ifort
  transfer_struct_c%bs_do_directbounce_only_f_byte_size = transfer_struct_c%bs_do_directbounce_only_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_wsabsa_output = .FALSE.
  bs_do_wsabsa_output_ptr => fortran_type_f%bs_do_wsabsa_output
  transfer_struct_c%bs_do_wsabsa_output = c_loc(bs_do_wsabsa_output_ptr)
  inquire(iolength=transfer_struct_c%bs_do_wsabsa_output_f_byte_size) fortran_type_f%bs_do_wsabsa_output
#ifdef ifort
  transfer_struct_c%bs_do_wsabsa_output_f_byte_size = transfer_struct_c%bs_do_wsabsa_output_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_wsa_scaling = .FALSE.
  bs_do_wsa_scaling_ptr => fortran_type_f%bs_do_wsa_scaling
  transfer_struct_c%bs_do_wsa_scaling = c_loc(bs_do_wsa_scaling_ptr)
  inquire(iolength=transfer_struct_c%bs_do_wsa_scaling_f_byte_size) fortran_type_f%bs_do_wsa_scaling
#ifdef ifort
  transfer_struct_c%bs_do_wsa_scaling_f_byte_size = transfer_struct_c%bs_do_wsa_scaling_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_bsa_scaling = .FALSE.
  bs_do_bsa_scaling_ptr => fortran_type_f%bs_do_bsa_scaling
  transfer_struct_c%bs_do_bsa_scaling = c_loc(bs_do_bsa_scaling_ptr)
  inquire(iolength=transfer_struct_c%bs_do_bsa_scaling_f_byte_size) fortran_type_f%bs_do_bsa_scaling
#ifdef ifort
  transfer_struct_c%bs_do_bsa_scaling_f_byte_size = transfer_struct_c%bs_do_bsa_scaling_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_wsa_value = 0_fpk
  bs_wsa_value_ptr => fortran_type_f%bs_wsa_value
  transfer_struct_c%bs_wsa_value = c_loc(bs_wsa_value_ptr)
  inquire(iolength=transfer_struct_c%bs_wsa_value_f_byte_size) fortran_type_f%bs_wsa_value
#ifdef ifort
  transfer_struct_c%bs_wsa_value_f_byte_size = transfer_struct_c%bs_wsa_value_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_bsa_value = 0_fpk
  bs_bsa_value_ptr => fortran_type_f%bs_bsa_value
  transfer_struct_c%bs_bsa_value = c_loc(bs_bsa_value_ptr)
  inquire(iolength=transfer_struct_c%bs_bsa_value_f_byte_size) fortran_type_f%bs_bsa_value
#ifdef ifort
  transfer_struct_c%bs_bsa_value_f_byte_size = transfer_struct_c%bs_bsa_value_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_newcmglint = .FALSE.
  bs_do_newcmglint_ptr => fortran_type_f%bs_do_newcmglint
  transfer_struct_c%bs_do_newcmglint = c_loc(bs_do_newcmglint_ptr)
  inquire(iolength=transfer_struct_c%bs_do_newcmglint_f_byte_size) fortran_type_f%bs_do_newcmglint
#ifdef ifort
  transfer_struct_c%bs_do_newcmglint_f_byte_size = transfer_struct_c%bs_do_newcmglint_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_salinity = 0_fpk
  bs_salinity_ptr => fortran_type_f%bs_salinity
  transfer_struct_c%bs_salinity = c_loc(bs_salinity_ptr)
  inquire(iolength=transfer_struct_c%bs_salinity_f_byte_size) fortran_type_f%bs_salinity
#ifdef ifort
  transfer_struct_c%bs_salinity_f_byte_size = transfer_struct_c%bs_salinity_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_wavelength = 0_fpk
  bs_wavelength_ptr => fortran_type_f%bs_wavelength
  transfer_struct_c%bs_wavelength = c_loc(bs_wavelength_ptr)
  inquire(iolength=transfer_struct_c%bs_wavelength_f_byte_size) fortran_type_f%bs_wavelength
#ifdef ifort
  transfer_struct_c%bs_wavelength_f_byte_size = transfer_struct_c%bs_wavelength_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_windspeed = 0_fpk
  bs_windspeed_ptr => fortran_type_f%bs_windspeed
  transfer_struct_c%bs_windspeed = c_loc(bs_windspeed_ptr)
  inquire(iolength=transfer_struct_c%bs_windspeed_f_byte_size) fortran_type_f%bs_windspeed
#ifdef ifort
  transfer_struct_c%bs_windspeed_f_byte_size = transfer_struct_c%bs_windspeed_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_winddir = 0_fpk
  bs_winddir_ptr => fortran_type_f%bs_winddir
  transfer_struct_c%bs_winddir = c_loc(bs_winddir_ptr(&
    lbound(fortran_type_f%bs_winddir,1)))
  inquire(iolength=transfer_struct_c%bs_winddir_f_byte_size) fortran_type_f%bs_winddir(&
    lbound(fortran_type_f%bs_winddir,1))
#ifdef ifort
  transfer_struct_c%bs_winddir_f_byte_size = transfer_struct_c%bs_winddir_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_winddir_f_shapes(1) = size(fortran_type_f%bs_winddir, 1)
  
  
  fortran_type_f%bs_do_glintshadow = .FALSE.
  bs_do_glintshadow_ptr => fortran_type_f%bs_do_glintshadow
  transfer_struct_c%bs_do_glintshadow = c_loc(bs_do_glintshadow_ptr)
  inquire(iolength=transfer_struct_c%bs_do_glintshadow_f_byte_size) fortran_type_f%bs_do_glintshadow
#ifdef ifort
  transfer_struct_c%bs_do_glintshadow_f_byte_size = transfer_struct_c%bs_do_glintshadow_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_foamoption = .FALSE.
  bs_do_foamoption_ptr => fortran_type_f%bs_do_foamoption
  transfer_struct_c%bs_do_foamoption = c_loc(bs_do_foamoption_ptr)
  inquire(iolength=transfer_struct_c%bs_do_foamoption_f_byte_size) fortran_type_f%bs_do_foamoption
#ifdef ifort
  transfer_struct_c%bs_do_foamoption_f_byte_size = transfer_struct_c%bs_do_foamoption_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_facetisotropy = .FALSE.
  bs_do_facetisotropy_ptr => fortran_type_f%bs_do_facetisotropy
  transfer_struct_c%bs_do_facetisotropy = c_loc(bs_do_facetisotropy_ptr)
  inquire(iolength=transfer_struct_c%bs_do_facetisotropy_f_byte_size) fortran_type_f%bs_do_facetisotropy
#ifdef ifort
  transfer_struct_c%bs_do_facetisotropy_f_byte_size = transfer_struct_c%bs_do_facetisotropy_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_glitter_msrcorr = .FALSE.
  bs_do_glitter_msrcorr_ptr => fortran_type_f%bs_do_glitter_msrcorr
  transfer_struct_c%bs_do_glitter_msrcorr = c_loc(bs_do_glitter_msrcorr_ptr)
  inquire(iolength=transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size) fortran_type_f%bs_do_glitter_msrcorr
#ifdef ifort
  transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size = transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_glitter_msrcorr_dbonly = .FALSE.
  bs_do_glitter_msrcorr_dbonly_ptr => fortran_type_f%bs_do_glitter_msrcorr_dbonly
  transfer_struct_c%bs_do_glitter_msrcorr_dbonly = c_loc(bs_do_glitter_msrcorr_dbonly_ptr)
  inquire(iolength=transfer_struct_c%bs_do_glitter_msrcorr_dbonly_f_byte_size) fortran_type_f%bs_do_glitter_msrcorr_dbonly
#ifdef ifort
  transfer_struct_c%bs_do_glitter_msrcorr_dbonly_f_byte_size = transfer_struct_c%bs_do_glitter_msrcorr_dbonly_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_glitter_msrcorr_order = 0
  bs_glitter_msrcorr_order_ptr => fortran_type_f%bs_glitter_msrcorr_order
  transfer_struct_c%bs_glitter_msrcorr_order = c_loc(bs_glitter_msrcorr_order_ptr)
  inquire(iolength=transfer_struct_c%bs_glitter_msrcorr_order_f_byte_size) fortran_type_f%bs_glitter_msrcorr_order
#ifdef ifort
  transfer_struct_c%bs_glitter_msrcorr_order_f_byte_size = transfer_struct_c%bs_glitter_msrcorr_order_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_glitter_msrcorr_nmuquad = 0
  bs_glitter_msrcorr_nmuquad_ptr => fortran_type_f%bs_glitter_msrcorr_nmuquad
  transfer_struct_c%bs_glitter_msrcorr_nmuquad = c_loc(bs_glitter_msrcorr_nmuquad_ptr)
  inquire(iolength=transfer_struct_c%bs_glitter_msrcorr_nmuquad_f_byte_size) fortran_type_f%bs_glitter_msrcorr_nmuquad
#ifdef ifort
  transfer_struct_c%bs_glitter_msrcorr_nmuquad_f_byte_size = transfer_struct_c%bs_glitter_msrcorr_nmuquad_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_glitter_msrcorr_nphiquad = 0
  bs_glitter_msrcorr_nphiquad_ptr => fortran_type_f%bs_glitter_msrcorr_nphiquad
  transfer_struct_c%bs_glitter_msrcorr_nphiquad = c_loc(bs_glitter_msrcorr_nphiquad_ptr)
  inquire(iolength=transfer_struct_c%bs_glitter_msrcorr_nphiquad_f_byte_size) fortran_type_f%bs_glitter_msrcorr_nphiquad
#ifdef ifort
  transfer_struct_c%bs_glitter_msrcorr_nphiquad_f_byte_size = transfer_struct_c%bs_glitter_msrcorr_nphiquad_f_byte_size * 4
#endif
  
  
  
end subroutine brdf_sup_inputs_c_init_only

subroutine brdf_sup_inputs_c_destroy(fortran_type_c) bind(C)
  use brdf_sup_inputs_def_m, only : brdf_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_sup_inputs_c_destroy

subroutine brdf_sup_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_inputs_def_m, only : brdf_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_sup_inputs), pointer :: fortran_type_f_from
  type(brdf_sup_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_do_brdf_surface = fortran_type_f_from%bs_do_brdf_surface
  fortran_type_f_to%bs_do_surface_emission = fortran_type_f_from%bs_do_surface_emission
  fortran_type_f_to%bs_do_solar_sources = fortran_type_f_from%bs_do_solar_sources
  fortran_type_f_to%bs_do_user_streams = fortran_type_f_from%bs_do_user_streams
  fortran_type_f_to%bs_do_user_obsgeoms = fortran_type_f_from%bs_do_user_obsgeoms
  fortran_type_f_to%bs_do_doublet_geometry = fortran_type_f_from%bs_do_doublet_geometry
  fortran_type_f_to%bs_nstreams = fortran_type_f_from%bs_nstreams
  fortran_type_f_to%bs_nbeams = fortran_type_f_from%bs_nbeams
  fortran_type_f_to%bs_beam_szas = fortran_type_f_from%bs_beam_szas
  fortran_type_f_to%bs_n_user_relazms = fortran_type_f_from%bs_n_user_relazms
  fortran_type_f_to%bs_user_relazms = fortran_type_f_from%bs_user_relazms
  fortran_type_f_to%bs_n_user_streams = fortran_type_f_from%bs_n_user_streams
  fortran_type_f_to%bs_user_angles_input = fortran_type_f_from%bs_user_angles_input
  fortran_type_f_to%bs_n_user_obsgeoms = fortran_type_f_from%bs_n_user_obsgeoms
  fortran_type_f_to%bs_user_obsgeoms = fortran_type_f_from%bs_user_obsgeoms
  fortran_type_f_to%bs_n_user_doublets = fortran_type_f_from%bs_n_user_doublets
  fortran_type_f_to%bs_user_doublets = fortran_type_f_from%bs_user_doublets
  fortran_type_f_to%bs_n_brdf_kernels = fortran_type_f_from%bs_n_brdf_kernels
  fortran_type_f_to%bs_brdf_names = fortran_type_f_from%bs_brdf_names
  fortran_type_f_to%bs_which_brdf = fortran_type_f_from%bs_which_brdf
  fortran_type_f_to%bs_n_brdf_parameters = fortran_type_f_from%bs_n_brdf_parameters
  fortran_type_f_to%bs_brdf_parameters = fortran_type_f_from%bs_brdf_parameters
  fortran_type_f_to%bs_lambertian_kernel_flag = fortran_type_f_from%bs_lambertian_kernel_flag
  fortran_type_f_to%bs_brdf_factors = fortran_type_f_from%bs_brdf_factors
  fortran_type_f_to%bs_nstreams_brdf = fortran_type_f_from%bs_nstreams_brdf
  fortran_type_f_to%bs_do_shadow_effect = fortran_type_f_from%bs_do_shadow_effect
  fortran_type_f_to%bs_do_directbounce_only = fortran_type_f_from%bs_do_directbounce_only
  fortran_type_f_to%bs_do_wsabsa_output = fortran_type_f_from%bs_do_wsabsa_output
  fortran_type_f_to%bs_do_wsa_scaling = fortran_type_f_from%bs_do_wsa_scaling
  fortran_type_f_to%bs_do_bsa_scaling = fortran_type_f_from%bs_do_bsa_scaling
  fortran_type_f_to%bs_wsa_value = fortran_type_f_from%bs_wsa_value
  fortran_type_f_to%bs_bsa_value = fortran_type_f_from%bs_bsa_value
  fortran_type_f_to%bs_do_newcmglint = fortran_type_f_from%bs_do_newcmglint
  fortran_type_f_to%bs_salinity = fortran_type_f_from%bs_salinity
  fortran_type_f_to%bs_wavelength = fortran_type_f_from%bs_wavelength
  fortran_type_f_to%bs_windspeed = fortran_type_f_from%bs_windspeed
  fortran_type_f_to%bs_winddir = fortran_type_f_from%bs_winddir
  fortran_type_f_to%bs_do_glintshadow = fortran_type_f_from%bs_do_glintshadow
  fortran_type_f_to%bs_do_foamoption = fortran_type_f_from%bs_do_foamoption
  fortran_type_f_to%bs_do_facetisotropy = fortran_type_f_from%bs_do_facetisotropy
  fortran_type_f_to%bs_do_glitter_msrcorr = fortran_type_f_from%bs_do_glitter_msrcorr
  fortran_type_f_to%bs_do_glitter_msrcorr_dbonly = fortran_type_f_from%bs_do_glitter_msrcorr_dbonly
  fortran_type_f_to%bs_glitter_msrcorr_order = fortran_type_f_from%bs_glitter_msrcorr_order
  fortran_type_f_to%bs_glitter_msrcorr_nmuquad = fortran_type_f_from%bs_glitter_msrcorr_nmuquad
  fortran_type_f_to%bs_glitter_msrcorr_nphiquad = fortran_type_f_from%bs_glitter_msrcorr_nphiquad
  

end subroutine brdf_sup_inputs_c_copy

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_sup_outputs

  type(brdf_sup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_sup_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_sup_outputs_c_alloc_init

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_sup_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m
  use lidort_pars_m

  type(brdf_sup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: bs_dbounce_brdfunc_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_user_brdf_f_ptr
  real(c_double), dimension(:), pointer :: bs_emissivity_ptr
  real(c_double), dimension(:), pointer :: bs_user_emissivity_ptr
  real(c_double), pointer :: bs_wsa_calculated_ptr
  real(c_double), dimension(:), pointer :: bs_wsa_kernels_ptr
  real(c_double), pointer :: bs_bsa_calculated_ptr
  real(c_double), dimension(:), pointer :: bs_bsa_kernels_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_dbounce_brdfunc = 0_fpk
  bs_dbounce_brdfunc_ptr => fortran_type_f%bs_dbounce_brdfunc
  transfer_struct_c%bs_dbounce_brdfunc = c_loc(bs_dbounce_brdfunc_ptr(&
    lbound(fortran_type_f%bs_dbounce_brdfunc,1),&
    lbound(fortran_type_f%bs_dbounce_brdfunc,2),&
    lbound(fortran_type_f%bs_dbounce_brdfunc,3)))
  inquire(iolength=transfer_struct_c%bs_dbounce_brdfunc_f_byte_size) fortran_type_f%bs_dbounce_brdfunc(&
    lbound(fortran_type_f%bs_dbounce_brdfunc,1),&
    lbound(fortran_type_f%bs_dbounce_brdfunc,2),&
    lbound(fortran_type_f%bs_dbounce_brdfunc,3))
#ifdef ifort
  transfer_struct_c%bs_dbounce_brdfunc_f_byte_size = transfer_struct_c%bs_dbounce_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_dbounce_brdfunc_f_shapes(1) = size(fortran_type_f%bs_dbounce_brdfunc, 1)
  transfer_struct_c%bs_dbounce_brdfunc_f_shapes(2) = size(fortran_type_f%bs_dbounce_brdfunc, 2)
  transfer_struct_c%bs_dbounce_brdfunc_f_shapes(3) = size(fortran_type_f%bs_dbounce_brdfunc, 3)
  
  
  fortran_type_f%bs_brdf_f_0 = 0_fpk
  bs_brdf_f_0_ptr => fortran_type_f%bs_brdf_f_0
  transfer_struct_c%bs_brdf_f_0 = c_loc(bs_brdf_f_0_ptr(&
    lbound(fortran_type_f%bs_brdf_f_0,1),&
    lbound(fortran_type_f%bs_brdf_f_0,2),&
    lbound(fortran_type_f%bs_brdf_f_0,3)))
  inquire(iolength=transfer_struct_c%bs_brdf_f_0_f_byte_size) fortran_type_f%bs_brdf_f_0(&
    lbound(fortran_type_f%bs_brdf_f_0,1),&
    lbound(fortran_type_f%bs_brdf_f_0,2),&
    lbound(fortran_type_f%bs_brdf_f_0,3))
#ifdef ifort
  transfer_struct_c%bs_brdf_f_0_f_byte_size = transfer_struct_c%bs_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_brdf_f_0_f_shapes(1) = size(fortran_type_f%bs_brdf_f_0, 1)
  transfer_struct_c%bs_brdf_f_0_f_shapes(2) = size(fortran_type_f%bs_brdf_f_0, 2)
  transfer_struct_c%bs_brdf_f_0_f_shapes(3) = size(fortran_type_f%bs_brdf_f_0, 3)
  
  
  fortran_type_f%bs_brdf_f = 0_fpk
  bs_brdf_f_ptr => fortran_type_f%bs_brdf_f
  transfer_struct_c%bs_brdf_f = c_loc(bs_brdf_f_ptr(&
    lbound(fortran_type_f%bs_brdf_f,1),&
    lbound(fortran_type_f%bs_brdf_f,2),&
    lbound(fortran_type_f%bs_brdf_f,3)))
  inquire(iolength=transfer_struct_c%bs_brdf_f_f_byte_size) fortran_type_f%bs_brdf_f(&
    lbound(fortran_type_f%bs_brdf_f,1),&
    lbound(fortran_type_f%bs_brdf_f,2),&
    lbound(fortran_type_f%bs_brdf_f,3))
#ifdef ifort
  transfer_struct_c%bs_brdf_f_f_byte_size = transfer_struct_c%bs_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_brdf_f_f_shapes(1) = size(fortran_type_f%bs_brdf_f, 1)
  transfer_struct_c%bs_brdf_f_f_shapes(2) = size(fortran_type_f%bs_brdf_f, 2)
  transfer_struct_c%bs_brdf_f_f_shapes(3) = size(fortran_type_f%bs_brdf_f, 3)
  
  
  fortran_type_f%bs_user_brdf_f_0 = 0_fpk
  bs_user_brdf_f_0_ptr => fortran_type_f%bs_user_brdf_f_0
  transfer_struct_c%bs_user_brdf_f_0 = c_loc(bs_user_brdf_f_0_ptr(&
    lbound(fortran_type_f%bs_user_brdf_f_0,1),&
    lbound(fortran_type_f%bs_user_brdf_f_0,2),&
    lbound(fortran_type_f%bs_user_brdf_f_0,3)))
  inquire(iolength=transfer_struct_c%bs_user_brdf_f_0_f_byte_size) fortran_type_f%bs_user_brdf_f_0(&
    lbound(fortran_type_f%bs_user_brdf_f_0,1),&
    lbound(fortran_type_f%bs_user_brdf_f_0,2),&
    lbound(fortran_type_f%bs_user_brdf_f_0,3))
#ifdef ifort
  transfer_struct_c%bs_user_brdf_f_0_f_byte_size = transfer_struct_c%bs_user_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_brdf_f_0_f_shapes(1) = size(fortran_type_f%bs_user_brdf_f_0, 1)
  transfer_struct_c%bs_user_brdf_f_0_f_shapes(2) = size(fortran_type_f%bs_user_brdf_f_0, 2)
  transfer_struct_c%bs_user_brdf_f_0_f_shapes(3) = size(fortran_type_f%bs_user_brdf_f_0, 3)
  
  
  fortran_type_f%bs_user_brdf_f = 0_fpk
  bs_user_brdf_f_ptr => fortran_type_f%bs_user_brdf_f
  transfer_struct_c%bs_user_brdf_f = c_loc(bs_user_brdf_f_ptr(&
    lbound(fortran_type_f%bs_user_brdf_f,1),&
    lbound(fortran_type_f%bs_user_brdf_f,2),&
    lbound(fortran_type_f%bs_user_brdf_f,3)))
  inquire(iolength=transfer_struct_c%bs_user_brdf_f_f_byte_size) fortran_type_f%bs_user_brdf_f(&
    lbound(fortran_type_f%bs_user_brdf_f,1),&
    lbound(fortran_type_f%bs_user_brdf_f,2),&
    lbound(fortran_type_f%bs_user_brdf_f,3))
#ifdef ifort
  transfer_struct_c%bs_user_brdf_f_f_byte_size = transfer_struct_c%bs_user_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_brdf_f_f_shapes(1) = size(fortran_type_f%bs_user_brdf_f, 1)
  transfer_struct_c%bs_user_brdf_f_f_shapes(2) = size(fortran_type_f%bs_user_brdf_f, 2)
  transfer_struct_c%bs_user_brdf_f_f_shapes(3) = size(fortran_type_f%bs_user_brdf_f, 3)
  
  
  fortran_type_f%bs_emissivity = 0_fpk
  bs_emissivity_ptr => fortran_type_f%bs_emissivity
  transfer_struct_c%bs_emissivity = c_loc(bs_emissivity_ptr(&
    lbound(fortran_type_f%bs_emissivity,1)))
  inquire(iolength=transfer_struct_c%bs_emissivity_f_byte_size) fortran_type_f%bs_emissivity(&
    lbound(fortran_type_f%bs_emissivity,1))
#ifdef ifort
  transfer_struct_c%bs_emissivity_f_byte_size = transfer_struct_c%bs_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_emissivity_f_shapes(1) = size(fortran_type_f%bs_emissivity, 1)
  
  
  fortran_type_f%bs_user_emissivity = 0_fpk
  bs_user_emissivity_ptr => fortran_type_f%bs_user_emissivity
  transfer_struct_c%bs_user_emissivity = c_loc(bs_user_emissivity_ptr(&
    lbound(fortran_type_f%bs_user_emissivity,1)))
  inquire(iolength=transfer_struct_c%bs_user_emissivity_f_byte_size) fortran_type_f%bs_user_emissivity(&
    lbound(fortran_type_f%bs_user_emissivity,1))
#ifdef ifort
  transfer_struct_c%bs_user_emissivity_f_byte_size = transfer_struct_c%bs_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_emissivity_f_shapes(1) = size(fortran_type_f%bs_user_emissivity, 1)
  
  
  fortran_type_f%bs_wsa_calculated = 0_fpk
  bs_wsa_calculated_ptr => fortran_type_f%bs_wsa_calculated
  transfer_struct_c%bs_wsa_calculated = c_loc(bs_wsa_calculated_ptr)
  inquire(iolength=transfer_struct_c%bs_wsa_calculated_f_byte_size) fortran_type_f%bs_wsa_calculated
#ifdef ifort
  transfer_struct_c%bs_wsa_calculated_f_byte_size = transfer_struct_c%bs_wsa_calculated_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_wsa_kernels = 0_fpk
  bs_wsa_kernels_ptr => fortran_type_f%bs_wsa_kernels
  transfer_struct_c%bs_wsa_kernels = c_loc(bs_wsa_kernels_ptr(&
    lbound(fortran_type_f%bs_wsa_kernels,1)))
  inquire(iolength=transfer_struct_c%bs_wsa_kernels_f_byte_size) fortran_type_f%bs_wsa_kernels(&
    lbound(fortran_type_f%bs_wsa_kernels,1))
#ifdef ifort
  transfer_struct_c%bs_wsa_kernels_f_byte_size = transfer_struct_c%bs_wsa_kernels_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_wsa_kernels_f_shapes(1) = size(fortran_type_f%bs_wsa_kernels, 1)
  
  
  fortran_type_f%bs_bsa_calculated = 0_fpk
  bs_bsa_calculated_ptr => fortran_type_f%bs_bsa_calculated
  transfer_struct_c%bs_bsa_calculated = c_loc(bs_bsa_calculated_ptr)
  inquire(iolength=transfer_struct_c%bs_bsa_calculated_f_byte_size) fortran_type_f%bs_bsa_calculated
#ifdef ifort
  transfer_struct_c%bs_bsa_calculated_f_byte_size = transfer_struct_c%bs_bsa_calculated_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_bsa_kernels = 0_fpk
  bs_bsa_kernels_ptr => fortran_type_f%bs_bsa_kernels
  transfer_struct_c%bs_bsa_kernels = c_loc(bs_bsa_kernels_ptr(&
    lbound(fortran_type_f%bs_bsa_kernels,1)))
  inquire(iolength=transfer_struct_c%bs_bsa_kernels_f_byte_size) fortran_type_f%bs_bsa_kernels(&
    lbound(fortran_type_f%bs_bsa_kernels,1))
#ifdef ifort
  transfer_struct_c%bs_bsa_kernels_f_byte_size = transfer_struct_c%bs_bsa_kernels_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_bsa_kernels_f_shapes(1) = size(fortran_type_f%bs_bsa_kernels, 1)
  
  
end subroutine brdf_sup_outputs_c_init_only

subroutine brdf_sup_outputs_c_destroy(fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_sup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_sup_outputs_c_destroy

subroutine brdf_sup_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_sup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_sup_outputs), pointer :: fortran_type_f_from
  type(brdf_sup_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_dbounce_brdfunc = fortran_type_f_from%bs_dbounce_brdfunc
  fortran_type_f_to%bs_brdf_f_0 = fortran_type_f_from%bs_brdf_f_0
  fortran_type_f_to%bs_brdf_f = fortran_type_f_from%bs_brdf_f
  fortran_type_f_to%bs_user_brdf_f_0 = fortran_type_f_from%bs_user_brdf_f_0
  fortran_type_f_to%bs_user_brdf_f = fortran_type_f_from%bs_user_brdf_f
  fortran_type_f_to%bs_emissivity = fortran_type_f_from%bs_emissivity
  fortran_type_f_to%bs_user_emissivity = fortran_type_f_from%bs_user_emissivity
  fortran_type_f_to%bs_wsa_calculated = fortran_type_f_from%bs_wsa_calculated
  fortran_type_f_to%bs_wsa_kernels = fortran_type_f_from%bs_wsa_kernels
  fortran_type_f_to%bs_bsa_calculated = fortran_type_f_from%bs_bsa_calculated
  fortran_type_f_to%bs_bsa_kernels = fortran_type_f_from%bs_bsa_kernels
  

end subroutine brdf_sup_outputs_c_copy

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_input_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_input_exception_handling

  type(brdf_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_input_exception_handling_c_alloc_init

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m
  use lidort_pars_m

  type(brdf_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  integer(c_int), pointer :: bs_status_inputread_ptr
  integer(c_int), pointer :: bs_ninputmessages_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_status_inputread = 0
  bs_status_inputread_ptr => fortran_type_f%bs_status_inputread
  transfer_struct_c%bs_status_inputread = c_loc(bs_status_inputread_ptr)
  inquire(iolength=transfer_struct_c%bs_status_inputread_f_byte_size) fortran_type_f%bs_status_inputread
#ifdef ifort
  transfer_struct_c%bs_status_inputread_f_byte_size = transfer_struct_c%bs_status_inputread_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_ninputmessages = 0
  bs_ninputmessages_ptr => fortran_type_f%bs_ninputmessages
  transfer_struct_c%bs_ninputmessages = c_loc(bs_ninputmessages_ptr)
  inquire(iolength=transfer_struct_c%bs_ninputmessages_f_byte_size) fortran_type_f%bs_ninputmessages
#ifdef ifort
  transfer_struct_c%bs_ninputmessages_f_byte_size = transfer_struct_c%bs_ninputmessages_f_byte_size * 4
#endif
  
  
  fortran_type_f%bs_inputmessages = ''
  transfer_struct_c%bs_inputmessages_f_len = len(fortran_type_f%bs_inputmessages)
  transfer_struct_c%bs_inputmessages_f_shapes(1) = size(fortran_type_f%bs_inputmessages, 1)
  
  fortran_type_f%bs_inputactions = ''
  transfer_struct_c%bs_inputactions_f_len = len(fortran_type_f%bs_inputactions)
  transfer_struct_c%bs_inputactions_f_shapes(1) = size(fortran_type_f%bs_inputactions, 1)
  
  
end subroutine brdf_input_exception_handling_c_init_only

subroutine brdf_input_exception_handling_c_destroy(fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_input_exception_handling_c_destroy

subroutine brdf_input_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_input_exception_handling), pointer :: fortran_type_f_from
  type(brdf_input_exception_handling), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_status_inputread = fortran_type_f_from%bs_status_inputread
  fortran_type_f_to%bs_ninputmessages = fortran_type_f_from%bs_ninputmessages
  fortran_type_f_to%bs_inputmessages = fortran_type_f_from%bs_inputmessages
  fortran_type_f_to%bs_inputactions = fortran_type_f_from%bs_inputactions
  

end subroutine brdf_input_exception_handling_c_copy

! Links to type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_output_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_output_exception_handling

  type(brdf_output_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_output_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_output_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_output_exception_handling_c_alloc_init

! Links to type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_output_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m
  use lidort_pars_m

  type(brdf_output_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_output_exception_handling), pointer :: fortran_type_f

  integer(c_int), pointer :: bs_status_output_ptr
  integer(c_int), pointer :: bs_noutputmessages_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_status_output = 0
  bs_status_output_ptr => fortran_type_f%bs_status_output
  transfer_struct_c%bs_status_output = c_loc(bs_status_output_ptr)
  inquire(iolength=transfer_struct_c%bs_status_output_f_byte_size) fortran_type_f%bs_status_output
#ifdef ifort
  transfer_struct_c%bs_status_output_f_byte_size = transfer_struct_c%bs_status_output_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_noutputmessages = 0
  bs_noutputmessages_ptr => fortran_type_f%bs_noutputmessages
  transfer_struct_c%bs_noutputmessages = c_loc(bs_noutputmessages_ptr)
  inquire(iolength=transfer_struct_c%bs_noutputmessages_f_byte_size) fortran_type_f%bs_noutputmessages
#ifdef ifort
  transfer_struct_c%bs_noutputmessages_f_byte_size = transfer_struct_c%bs_noutputmessages_f_byte_size * 4
#endif
  
  
  fortran_type_f%bs_outputmessages = ''
  transfer_struct_c%bs_outputmessages_f_len = len(fortran_type_f%bs_outputmessages)
  transfer_struct_c%bs_outputmessages_f_shapes(1) = size(fortran_type_f%bs_outputmessages, 1)
  
  
end subroutine brdf_output_exception_handling_c_init_only

subroutine brdf_output_exception_handling_c_destroy(fortran_type_c) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_output_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_output_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_output_exception_handling_c_destroy

subroutine brdf_output_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_output_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_output_exception_handling), pointer :: fortran_type_f_from
  type(brdf_output_exception_handling), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_status_output = fortran_type_f_from%bs_status_output
  fortran_type_f_to%bs_noutputmessages = fortran_type_f_from%bs_noutputmessages
  fortran_type_f_to%bs_outputmessages = fortran_type_f_from%bs_outputmessages
  

end subroutine brdf_output_exception_handling_c_copy

! Links to type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.F90"
! Allocs and initializes type
subroutine sleave_sup_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use sleave_sup_inputs_def_m, only : sleave_sup_inputs

  type(sleave_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(sleave_sup_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call sleave_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine sleave_sup_inputs_c_alloc_init

! Links to type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.F90"
! Initializes only with no allocation
subroutine sleave_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use sleave_sup_inputs_def_m
  use lidort_pars_m

  type(sleave_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(sleave_sup_inputs), pointer :: fortran_type_f

  logical(kind=4), pointer :: sl_do_sleaving_ptr
  logical(kind=4), pointer :: sl_do_isotropic_ptr
  logical(kind=4), pointer :: sl_do_roughsurface_ptr
  logical(kind=4), pointer :: sl_do_exact_ptr
  logical(kind=4), pointer :: sl_do_exactonly_ptr
  logical(kind=4), pointer :: sl_do_fluorescence_ptr
  logical(kind=4), pointer :: sl_do_solar_sources_ptr
  logical(kind=4), pointer :: sl_do_user_streams_ptr
  logical(kind=4), pointer :: sl_do_user_obsgeoms_ptr
  logical(kind=4), pointer :: sl_do_doublet_geometry_ptr
  integer(c_int), pointer :: sl_nstreams_ptr
  integer(c_int), pointer :: sl_nbeams_ptr
  real(c_double), dimension(:), pointer :: sl_beam_szas_ptr
  integer(c_int), pointer :: sl_n_user_relazms_ptr
  real(c_double), dimension(:), pointer :: sl_user_relazms_ptr
  integer(c_int), pointer :: sl_n_user_streams_ptr
  real(c_double), dimension(:), pointer :: sl_user_angles_input_ptr
  integer(c_int), pointer :: sl_n_user_obsgeoms_ptr
  real(c_double), dimension(:,:), pointer :: sl_user_obsgeoms_ptr
  integer(c_int), pointer :: sl_n_user_doublets_ptr
  real(c_double), dimension(:,:), pointer :: sl_user_doublets_ptr
  real(c_double), pointer :: sl_salinity_ptr
  real(c_double), pointer :: sl_chlorconc_ptr
  real(c_double), pointer :: sl_wavelength_ptr
  logical(kind=4), pointer :: sl_azimuthdep_ptr
  logical(kind=4), pointer :: sl_do_fourier_output_ptr
  real(c_double), pointer :: sl_windspeed_ptr
  real(c_double), dimension(:), pointer :: sl_winddir_ptr
  logical(kind=4), pointer :: sl_do_glintshadow_ptr
  logical(kind=4), pointer :: sl_do_foamoption_ptr
  logical(kind=4), pointer :: sl_do_facetisotropy_ptr
  real(c_double), pointer :: sl_fl_wavelength_ptr
  real(c_double), pointer :: sl_fl_latitude_ptr
  real(c_double), pointer :: sl_fl_longitude_ptr
  integer(c_int), dimension(:), pointer :: sl_fl_epoch_ptr
  real(c_double), pointer :: sl_fl_amplitude755_ptr
  logical(kind=4), pointer :: sl_fl_do_datagaussian_ptr
  real(c_double), dimension(:,:), pointer :: sl_fl_inputgaussians_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%sl_do_sleaving = .FALSE.
  sl_do_sleaving_ptr => fortran_type_f%sl_do_sleaving
  transfer_struct_c%sl_do_sleaving = c_loc(sl_do_sleaving_ptr)
  inquire(iolength=transfer_struct_c%sl_do_sleaving_f_byte_size) fortran_type_f%sl_do_sleaving
#ifdef ifort
  transfer_struct_c%sl_do_sleaving_f_byte_size = transfer_struct_c%sl_do_sleaving_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_isotropic = .FALSE.
  sl_do_isotropic_ptr => fortran_type_f%sl_do_isotropic
  transfer_struct_c%sl_do_isotropic = c_loc(sl_do_isotropic_ptr)
  inquire(iolength=transfer_struct_c%sl_do_isotropic_f_byte_size) fortran_type_f%sl_do_isotropic
#ifdef ifort
  transfer_struct_c%sl_do_isotropic_f_byte_size = transfer_struct_c%sl_do_isotropic_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_roughsurface = .FALSE.
  sl_do_roughsurface_ptr => fortran_type_f%sl_do_roughsurface
  transfer_struct_c%sl_do_roughsurface = c_loc(sl_do_roughsurface_ptr)
  inquire(iolength=transfer_struct_c%sl_do_roughsurface_f_byte_size) fortran_type_f%sl_do_roughsurface
#ifdef ifort
  transfer_struct_c%sl_do_roughsurface_f_byte_size = transfer_struct_c%sl_do_roughsurface_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_exact = .FALSE.
  sl_do_exact_ptr => fortran_type_f%sl_do_exact
  transfer_struct_c%sl_do_exact = c_loc(sl_do_exact_ptr)
  inquire(iolength=transfer_struct_c%sl_do_exact_f_byte_size) fortran_type_f%sl_do_exact
#ifdef ifort
  transfer_struct_c%sl_do_exact_f_byte_size = transfer_struct_c%sl_do_exact_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_exactonly = .FALSE.
  sl_do_exactonly_ptr => fortran_type_f%sl_do_exactonly
  transfer_struct_c%sl_do_exactonly = c_loc(sl_do_exactonly_ptr)
  inquire(iolength=transfer_struct_c%sl_do_exactonly_f_byte_size) fortran_type_f%sl_do_exactonly
#ifdef ifort
  transfer_struct_c%sl_do_exactonly_f_byte_size = transfer_struct_c%sl_do_exactonly_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_fluorescence = .FALSE.
  sl_do_fluorescence_ptr => fortran_type_f%sl_do_fluorescence
  transfer_struct_c%sl_do_fluorescence = c_loc(sl_do_fluorescence_ptr)
  inquire(iolength=transfer_struct_c%sl_do_fluorescence_f_byte_size) fortran_type_f%sl_do_fluorescence
#ifdef ifort
  transfer_struct_c%sl_do_fluorescence_f_byte_size = transfer_struct_c%sl_do_fluorescence_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_solar_sources = .FALSE.
  sl_do_solar_sources_ptr => fortran_type_f%sl_do_solar_sources
  transfer_struct_c%sl_do_solar_sources = c_loc(sl_do_solar_sources_ptr)
  inquire(iolength=transfer_struct_c%sl_do_solar_sources_f_byte_size) fortran_type_f%sl_do_solar_sources
#ifdef ifort
  transfer_struct_c%sl_do_solar_sources_f_byte_size = transfer_struct_c%sl_do_solar_sources_f_byte_size * 4
#endif
  
  
  fortran_type_f%sl_sleave_datapath = ''
  transfer_struct_c%sl_sleave_datapath_f_len = len(fortran_type_f%sl_sleave_datapath)
  
  
  fortran_type_f%sl_do_user_streams = .FALSE.
  sl_do_user_streams_ptr => fortran_type_f%sl_do_user_streams
  transfer_struct_c%sl_do_user_streams = c_loc(sl_do_user_streams_ptr)
  inquire(iolength=transfer_struct_c%sl_do_user_streams_f_byte_size) fortran_type_f%sl_do_user_streams
#ifdef ifort
  transfer_struct_c%sl_do_user_streams_f_byte_size = transfer_struct_c%sl_do_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_user_obsgeoms = .FALSE.
  sl_do_user_obsgeoms_ptr => fortran_type_f%sl_do_user_obsgeoms
  transfer_struct_c%sl_do_user_obsgeoms = c_loc(sl_do_user_obsgeoms_ptr)
  inquire(iolength=transfer_struct_c%sl_do_user_obsgeoms_f_byte_size) fortran_type_f%sl_do_user_obsgeoms
#ifdef ifort
  transfer_struct_c%sl_do_user_obsgeoms_f_byte_size = transfer_struct_c%sl_do_user_obsgeoms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_doublet_geometry = .FALSE.
  sl_do_doublet_geometry_ptr => fortran_type_f%sl_do_doublet_geometry
  transfer_struct_c%sl_do_doublet_geometry = c_loc(sl_do_doublet_geometry_ptr)
  inquire(iolength=transfer_struct_c%sl_do_doublet_geometry_f_byte_size) fortran_type_f%sl_do_doublet_geometry
#ifdef ifort
  transfer_struct_c%sl_do_doublet_geometry_f_byte_size = transfer_struct_c%sl_do_doublet_geometry_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_nstreams = 0
  sl_nstreams_ptr => fortran_type_f%sl_nstreams
  transfer_struct_c%sl_nstreams = c_loc(sl_nstreams_ptr)
  inquire(iolength=transfer_struct_c%sl_nstreams_f_byte_size) fortran_type_f%sl_nstreams
#ifdef ifort
  transfer_struct_c%sl_nstreams_f_byte_size = transfer_struct_c%sl_nstreams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_nbeams = 0
  sl_nbeams_ptr => fortran_type_f%sl_nbeams
  transfer_struct_c%sl_nbeams = c_loc(sl_nbeams_ptr)
  inquire(iolength=transfer_struct_c%sl_nbeams_f_byte_size) fortran_type_f%sl_nbeams
#ifdef ifort
  transfer_struct_c%sl_nbeams_f_byte_size = transfer_struct_c%sl_nbeams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_beam_szas = 0_fpk
  sl_beam_szas_ptr => fortran_type_f%sl_beam_szas
  transfer_struct_c%sl_beam_szas = c_loc(sl_beam_szas_ptr(&
    lbound(fortran_type_f%sl_beam_szas,1)))
  inquire(iolength=transfer_struct_c%sl_beam_szas_f_byte_size) fortran_type_f%sl_beam_szas(&
    lbound(fortran_type_f%sl_beam_szas,1))
#ifdef ifort
  transfer_struct_c%sl_beam_szas_f_byte_size = transfer_struct_c%sl_beam_szas_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_beam_szas_f_shapes(1) = size(fortran_type_f%sl_beam_szas, 1)
  
  
  fortran_type_f%sl_n_user_relazms = 0
  sl_n_user_relazms_ptr => fortran_type_f%sl_n_user_relazms
  transfer_struct_c%sl_n_user_relazms = c_loc(sl_n_user_relazms_ptr)
  inquire(iolength=transfer_struct_c%sl_n_user_relazms_f_byte_size) fortran_type_f%sl_n_user_relazms
#ifdef ifort
  transfer_struct_c%sl_n_user_relazms_f_byte_size = transfer_struct_c%sl_n_user_relazms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_user_relazms = 0_fpk
  sl_user_relazms_ptr => fortran_type_f%sl_user_relazms
  transfer_struct_c%sl_user_relazms = c_loc(sl_user_relazms_ptr(&
    lbound(fortran_type_f%sl_user_relazms,1)))
  inquire(iolength=transfer_struct_c%sl_user_relazms_f_byte_size) fortran_type_f%sl_user_relazms(&
    lbound(fortran_type_f%sl_user_relazms,1))
#ifdef ifort
  transfer_struct_c%sl_user_relazms_f_byte_size = transfer_struct_c%sl_user_relazms_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_user_relazms_f_shapes(1) = size(fortran_type_f%sl_user_relazms, 1)
  
  
  fortran_type_f%sl_n_user_streams = 0
  sl_n_user_streams_ptr => fortran_type_f%sl_n_user_streams
  transfer_struct_c%sl_n_user_streams = c_loc(sl_n_user_streams_ptr)
  inquire(iolength=transfer_struct_c%sl_n_user_streams_f_byte_size) fortran_type_f%sl_n_user_streams
#ifdef ifort
  transfer_struct_c%sl_n_user_streams_f_byte_size = transfer_struct_c%sl_n_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_user_angles_input = 0_fpk
  sl_user_angles_input_ptr => fortran_type_f%sl_user_angles_input
  transfer_struct_c%sl_user_angles_input = c_loc(sl_user_angles_input_ptr(&
    lbound(fortran_type_f%sl_user_angles_input,1)))
  inquire(iolength=transfer_struct_c%sl_user_angles_input_f_byte_size) fortran_type_f%sl_user_angles_input(&
    lbound(fortran_type_f%sl_user_angles_input,1))
#ifdef ifort
  transfer_struct_c%sl_user_angles_input_f_byte_size = transfer_struct_c%sl_user_angles_input_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_user_angles_input_f_shapes(1) = size(fortran_type_f%sl_user_angles_input, 1)
  
  
  fortran_type_f%sl_n_user_obsgeoms = 0
  sl_n_user_obsgeoms_ptr => fortran_type_f%sl_n_user_obsgeoms
  transfer_struct_c%sl_n_user_obsgeoms = c_loc(sl_n_user_obsgeoms_ptr)
  inquire(iolength=transfer_struct_c%sl_n_user_obsgeoms_f_byte_size) fortran_type_f%sl_n_user_obsgeoms
#ifdef ifort
  transfer_struct_c%sl_n_user_obsgeoms_f_byte_size = transfer_struct_c%sl_n_user_obsgeoms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_user_obsgeoms = 0_fpk
  sl_user_obsgeoms_ptr => fortran_type_f%sl_user_obsgeoms
  transfer_struct_c%sl_user_obsgeoms = c_loc(sl_user_obsgeoms_ptr(&
    lbound(fortran_type_f%sl_user_obsgeoms,1),&
    lbound(fortran_type_f%sl_user_obsgeoms,2)))
  inquire(iolength=transfer_struct_c%sl_user_obsgeoms_f_byte_size) fortran_type_f%sl_user_obsgeoms(&
    lbound(fortran_type_f%sl_user_obsgeoms,1),&
    lbound(fortran_type_f%sl_user_obsgeoms,2))
#ifdef ifort
  transfer_struct_c%sl_user_obsgeoms_f_byte_size = transfer_struct_c%sl_user_obsgeoms_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_user_obsgeoms_f_shapes(1) = size(fortran_type_f%sl_user_obsgeoms, 1)
  transfer_struct_c%sl_user_obsgeoms_f_shapes(2) = size(fortran_type_f%sl_user_obsgeoms, 2)
  
  
  fortran_type_f%sl_n_user_doublets = 0
  sl_n_user_doublets_ptr => fortran_type_f%sl_n_user_doublets
  transfer_struct_c%sl_n_user_doublets = c_loc(sl_n_user_doublets_ptr)
  inquire(iolength=transfer_struct_c%sl_n_user_doublets_f_byte_size) fortran_type_f%sl_n_user_doublets
#ifdef ifort
  transfer_struct_c%sl_n_user_doublets_f_byte_size = transfer_struct_c%sl_n_user_doublets_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_user_doublets = 0_fpk
  sl_user_doublets_ptr => fortran_type_f%sl_user_doublets
  transfer_struct_c%sl_user_doublets = c_loc(sl_user_doublets_ptr(&
    lbound(fortran_type_f%sl_user_doublets,1),&
    lbound(fortran_type_f%sl_user_doublets,2)))
  inquire(iolength=transfer_struct_c%sl_user_doublets_f_byte_size) fortran_type_f%sl_user_doublets(&
    lbound(fortran_type_f%sl_user_doublets,1),&
    lbound(fortran_type_f%sl_user_doublets,2))
#ifdef ifort
  transfer_struct_c%sl_user_doublets_f_byte_size = transfer_struct_c%sl_user_doublets_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_user_doublets_f_shapes(1) = size(fortran_type_f%sl_user_doublets, 1)
  transfer_struct_c%sl_user_doublets_f_shapes(2) = size(fortran_type_f%sl_user_doublets, 2)
  
  
  fortran_type_f%sl_salinity = 0_fpk
  sl_salinity_ptr => fortran_type_f%sl_salinity
  transfer_struct_c%sl_salinity = c_loc(sl_salinity_ptr)
  inquire(iolength=transfer_struct_c%sl_salinity_f_byte_size) fortran_type_f%sl_salinity
#ifdef ifort
  transfer_struct_c%sl_salinity_f_byte_size = transfer_struct_c%sl_salinity_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_chlorconc = 0_fpk
  sl_chlorconc_ptr => fortran_type_f%sl_chlorconc
  transfer_struct_c%sl_chlorconc = c_loc(sl_chlorconc_ptr)
  inquire(iolength=transfer_struct_c%sl_chlorconc_f_byte_size) fortran_type_f%sl_chlorconc
#ifdef ifort
  transfer_struct_c%sl_chlorconc_f_byte_size = transfer_struct_c%sl_chlorconc_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_wavelength = 0_fpk
  sl_wavelength_ptr => fortran_type_f%sl_wavelength
  transfer_struct_c%sl_wavelength = c_loc(sl_wavelength_ptr)
  inquire(iolength=transfer_struct_c%sl_wavelength_f_byte_size) fortran_type_f%sl_wavelength
#ifdef ifort
  transfer_struct_c%sl_wavelength_f_byte_size = transfer_struct_c%sl_wavelength_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_azimuthdep = .FALSE.
  sl_azimuthdep_ptr => fortran_type_f%sl_azimuthdep
  transfer_struct_c%sl_azimuthdep = c_loc(sl_azimuthdep_ptr)
  inquire(iolength=transfer_struct_c%sl_azimuthdep_f_byte_size) fortran_type_f%sl_azimuthdep
#ifdef ifort
  transfer_struct_c%sl_azimuthdep_f_byte_size = transfer_struct_c%sl_azimuthdep_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_fourier_output = .FALSE.
  sl_do_fourier_output_ptr => fortran_type_f%sl_do_fourier_output
  transfer_struct_c%sl_do_fourier_output = c_loc(sl_do_fourier_output_ptr)
  inquire(iolength=transfer_struct_c%sl_do_fourier_output_f_byte_size) fortran_type_f%sl_do_fourier_output
#ifdef ifort
  transfer_struct_c%sl_do_fourier_output_f_byte_size = transfer_struct_c%sl_do_fourier_output_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_windspeed = 0_fpk
  sl_windspeed_ptr => fortran_type_f%sl_windspeed
  transfer_struct_c%sl_windspeed = c_loc(sl_windspeed_ptr)
  inquire(iolength=transfer_struct_c%sl_windspeed_f_byte_size) fortran_type_f%sl_windspeed
#ifdef ifort
  transfer_struct_c%sl_windspeed_f_byte_size = transfer_struct_c%sl_windspeed_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_winddir = 0_fpk
  sl_winddir_ptr => fortran_type_f%sl_winddir
  transfer_struct_c%sl_winddir = c_loc(sl_winddir_ptr(&
    lbound(fortran_type_f%sl_winddir,1)))
  inquire(iolength=transfer_struct_c%sl_winddir_f_byte_size) fortran_type_f%sl_winddir(&
    lbound(fortran_type_f%sl_winddir,1))
#ifdef ifort
  transfer_struct_c%sl_winddir_f_byte_size = transfer_struct_c%sl_winddir_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_winddir_f_shapes(1) = size(fortran_type_f%sl_winddir, 1)
  
  
  fortran_type_f%sl_do_glintshadow = .FALSE.
  sl_do_glintshadow_ptr => fortran_type_f%sl_do_glintshadow
  transfer_struct_c%sl_do_glintshadow = c_loc(sl_do_glintshadow_ptr)
  inquire(iolength=transfer_struct_c%sl_do_glintshadow_f_byte_size) fortran_type_f%sl_do_glintshadow
#ifdef ifort
  transfer_struct_c%sl_do_glintshadow_f_byte_size = transfer_struct_c%sl_do_glintshadow_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_foamoption = .FALSE.
  sl_do_foamoption_ptr => fortran_type_f%sl_do_foamoption
  transfer_struct_c%sl_do_foamoption = c_loc(sl_do_foamoption_ptr)
  inquire(iolength=transfer_struct_c%sl_do_foamoption_f_byte_size) fortran_type_f%sl_do_foamoption
#ifdef ifort
  transfer_struct_c%sl_do_foamoption_f_byte_size = transfer_struct_c%sl_do_foamoption_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_do_facetisotropy = .FALSE.
  sl_do_facetisotropy_ptr => fortran_type_f%sl_do_facetisotropy
  transfer_struct_c%sl_do_facetisotropy = c_loc(sl_do_facetisotropy_ptr)
  inquire(iolength=transfer_struct_c%sl_do_facetisotropy_f_byte_size) fortran_type_f%sl_do_facetisotropy
#ifdef ifort
  transfer_struct_c%sl_do_facetisotropy_f_byte_size = transfer_struct_c%sl_do_facetisotropy_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_wavelength = 0_fpk
  sl_fl_wavelength_ptr => fortran_type_f%sl_fl_wavelength
  transfer_struct_c%sl_fl_wavelength = c_loc(sl_fl_wavelength_ptr)
  inquire(iolength=transfer_struct_c%sl_fl_wavelength_f_byte_size) fortran_type_f%sl_fl_wavelength
#ifdef ifort
  transfer_struct_c%sl_fl_wavelength_f_byte_size = transfer_struct_c%sl_fl_wavelength_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_latitude = 0_fpk
  sl_fl_latitude_ptr => fortran_type_f%sl_fl_latitude
  transfer_struct_c%sl_fl_latitude = c_loc(sl_fl_latitude_ptr)
  inquire(iolength=transfer_struct_c%sl_fl_latitude_f_byte_size) fortran_type_f%sl_fl_latitude
#ifdef ifort
  transfer_struct_c%sl_fl_latitude_f_byte_size = transfer_struct_c%sl_fl_latitude_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_longitude = 0_fpk
  sl_fl_longitude_ptr => fortran_type_f%sl_fl_longitude
  transfer_struct_c%sl_fl_longitude = c_loc(sl_fl_longitude_ptr)
  inquire(iolength=transfer_struct_c%sl_fl_longitude_f_byte_size) fortran_type_f%sl_fl_longitude
#ifdef ifort
  transfer_struct_c%sl_fl_longitude_f_byte_size = transfer_struct_c%sl_fl_longitude_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_epoch = 0
  sl_fl_epoch_ptr => fortran_type_f%sl_fl_epoch
  transfer_struct_c%sl_fl_epoch = c_loc(sl_fl_epoch_ptr(&
    lbound(fortran_type_f%sl_fl_epoch,1)))
  inquire(iolength=transfer_struct_c%sl_fl_epoch_f_byte_size) fortran_type_f%sl_fl_epoch(&
    lbound(fortran_type_f%sl_fl_epoch,1))
#ifdef ifort
  transfer_struct_c%sl_fl_epoch_f_byte_size = transfer_struct_c%sl_fl_epoch_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_fl_epoch_f_shapes(1) = size(fortran_type_f%sl_fl_epoch, 1)
  
  
  fortran_type_f%sl_fl_amplitude755 = 0_fpk
  sl_fl_amplitude755_ptr => fortran_type_f%sl_fl_amplitude755
  transfer_struct_c%sl_fl_amplitude755 = c_loc(sl_fl_amplitude755_ptr)
  inquire(iolength=transfer_struct_c%sl_fl_amplitude755_f_byte_size) fortran_type_f%sl_fl_amplitude755
#ifdef ifort
  transfer_struct_c%sl_fl_amplitude755_f_byte_size = transfer_struct_c%sl_fl_amplitude755_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_do_datagaussian = .FALSE.
  sl_fl_do_datagaussian_ptr => fortran_type_f%sl_fl_do_datagaussian
  transfer_struct_c%sl_fl_do_datagaussian = c_loc(sl_fl_do_datagaussian_ptr)
  inquire(iolength=transfer_struct_c%sl_fl_do_datagaussian_f_byte_size) fortran_type_f%sl_fl_do_datagaussian
#ifdef ifort
  transfer_struct_c%sl_fl_do_datagaussian_f_byte_size = transfer_struct_c%sl_fl_do_datagaussian_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%sl_fl_inputgaussians = 0_fpk
  sl_fl_inputgaussians_ptr => fortran_type_f%sl_fl_inputgaussians
  transfer_struct_c%sl_fl_inputgaussians = c_loc(sl_fl_inputgaussians_ptr(&
    lbound(fortran_type_f%sl_fl_inputgaussians,1),&
    lbound(fortran_type_f%sl_fl_inputgaussians,2)))
  inquire(iolength=transfer_struct_c%sl_fl_inputgaussians_f_byte_size) fortran_type_f%sl_fl_inputgaussians(&
    lbound(fortran_type_f%sl_fl_inputgaussians,1),&
    lbound(fortran_type_f%sl_fl_inputgaussians,2))
#ifdef ifort
  transfer_struct_c%sl_fl_inputgaussians_f_byte_size = transfer_struct_c%sl_fl_inputgaussians_f_byte_size * 4
#endif
  
  transfer_struct_c%sl_fl_inputgaussians_f_shapes(1) = size(fortran_type_f%sl_fl_inputgaussians, 1)
  transfer_struct_c%sl_fl_inputgaussians_f_shapes(2) = size(fortran_type_f%sl_fl_inputgaussians, 2)
  
  
end subroutine sleave_sup_inputs_c_init_only

subroutine sleave_sup_inputs_c_destroy(fortran_type_c) bind(C)
  use sleave_sup_inputs_def_m, only : sleave_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(sleave_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine sleave_sup_inputs_c_destroy

subroutine sleave_sup_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use sleave_sup_inputs_def_m, only : sleave_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(sleave_sup_inputs), pointer :: fortran_type_f_from
  type(sleave_sup_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%sl_do_sleaving = fortran_type_f_from%sl_do_sleaving
  fortran_type_f_to%sl_do_isotropic = fortran_type_f_from%sl_do_isotropic
  fortran_type_f_to%sl_do_roughsurface = fortran_type_f_from%sl_do_roughsurface
  fortran_type_f_to%sl_do_exact = fortran_type_f_from%sl_do_exact
  fortran_type_f_to%sl_do_exactonly = fortran_type_f_from%sl_do_exactonly
  fortran_type_f_to%sl_do_fluorescence = fortran_type_f_from%sl_do_fluorescence
  fortran_type_f_to%sl_do_solar_sources = fortran_type_f_from%sl_do_solar_sources
  fortran_type_f_to%sl_sleave_datapath = fortran_type_f_from%sl_sleave_datapath
  fortran_type_f_to%sl_do_user_streams = fortran_type_f_from%sl_do_user_streams
  fortran_type_f_to%sl_do_user_obsgeoms = fortran_type_f_from%sl_do_user_obsgeoms
  fortran_type_f_to%sl_do_doublet_geometry = fortran_type_f_from%sl_do_doublet_geometry
  fortran_type_f_to%sl_nstreams = fortran_type_f_from%sl_nstreams
  fortran_type_f_to%sl_nbeams = fortran_type_f_from%sl_nbeams
  fortran_type_f_to%sl_beam_szas = fortran_type_f_from%sl_beam_szas
  fortran_type_f_to%sl_n_user_relazms = fortran_type_f_from%sl_n_user_relazms
  fortran_type_f_to%sl_user_relazms = fortran_type_f_from%sl_user_relazms
  fortran_type_f_to%sl_n_user_streams = fortran_type_f_from%sl_n_user_streams
  fortran_type_f_to%sl_user_angles_input = fortran_type_f_from%sl_user_angles_input
  fortran_type_f_to%sl_n_user_obsgeoms = fortran_type_f_from%sl_n_user_obsgeoms
  fortran_type_f_to%sl_user_obsgeoms = fortran_type_f_from%sl_user_obsgeoms
  fortran_type_f_to%sl_n_user_doublets = fortran_type_f_from%sl_n_user_doublets
  fortran_type_f_to%sl_user_doublets = fortran_type_f_from%sl_user_doublets
  fortran_type_f_to%sl_salinity = fortran_type_f_from%sl_salinity
  fortran_type_f_to%sl_chlorconc = fortran_type_f_from%sl_chlorconc
  fortran_type_f_to%sl_wavelength = fortran_type_f_from%sl_wavelength
  fortran_type_f_to%sl_azimuthdep = fortran_type_f_from%sl_azimuthdep
  fortran_type_f_to%sl_do_fourier_output = fortran_type_f_from%sl_do_fourier_output
  fortran_type_f_to%sl_windspeed = fortran_type_f_from%sl_windspeed
  fortran_type_f_to%sl_winddir = fortran_type_f_from%sl_winddir
  fortran_type_f_to%sl_do_glintshadow = fortran_type_f_from%sl_do_glintshadow
  fortran_type_f_to%sl_do_foamoption = fortran_type_f_from%sl_do_foamoption
  fortran_type_f_to%sl_do_facetisotropy = fortran_type_f_from%sl_do_facetisotropy
  fortran_type_f_to%sl_fl_wavelength = fortran_type_f_from%sl_fl_wavelength
  fortran_type_f_to%sl_fl_latitude = fortran_type_f_from%sl_fl_latitude
  fortran_type_f_to%sl_fl_longitude = fortran_type_f_from%sl_fl_longitude
  fortran_type_f_to%sl_fl_epoch = fortran_type_f_from%sl_fl_epoch
  fortran_type_f_to%sl_fl_amplitude755 = fortran_type_f_from%sl_fl_amplitude755
  fortran_type_f_to%sl_fl_do_datagaussian = fortran_type_f_from%sl_fl_do_datagaussian
  fortran_type_f_to%sl_fl_inputgaussians = fortran_type_f_from%sl_fl_inputgaussians
  

end subroutine sleave_sup_inputs_c_copy

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lincontrol_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lincontrol

  type(lidort_fixed_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_lincontrol_c_alloc_init

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  logical(kind=4), dimension(:), pointer :: ts_layer_vary_flag_ptr
  integer(c_int), dimension(:), pointer :: ts_layer_vary_number_ptr
  integer(c_int), pointer :: ts_n_totalcolumn_wfs_ptr
  integer(c_int), pointer :: ts_n_surface_wfs_ptr
  integer(c_int), pointer :: ts_n_sleave_wfs_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_layer_vary_flag = .FALSE.
  ts_layer_vary_flag_ptr => fortran_type_f%ts_layer_vary_flag
  transfer_struct_c%ts_layer_vary_flag = c_loc(ts_layer_vary_flag_ptr(&
    lbound(fortran_type_f%ts_layer_vary_flag,1)))
  inquire(iolength=transfer_struct_c%ts_layer_vary_flag_f_byte_size) fortran_type_f%ts_layer_vary_flag(&
    lbound(fortran_type_f%ts_layer_vary_flag,1))
#ifdef ifort
  transfer_struct_c%ts_layer_vary_flag_f_byte_size = transfer_struct_c%ts_layer_vary_flag_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_layer_vary_flag_f_shapes(1) = size(fortran_type_f%ts_layer_vary_flag, 1)
  
  
  fortran_type_f%ts_layer_vary_number = 0
  ts_layer_vary_number_ptr => fortran_type_f%ts_layer_vary_number
  transfer_struct_c%ts_layer_vary_number = c_loc(ts_layer_vary_number_ptr(&
    lbound(fortran_type_f%ts_layer_vary_number,1)))
  inquire(iolength=transfer_struct_c%ts_layer_vary_number_f_byte_size) fortran_type_f%ts_layer_vary_number(&
    lbound(fortran_type_f%ts_layer_vary_number,1))
#ifdef ifort
  transfer_struct_c%ts_layer_vary_number_f_byte_size = transfer_struct_c%ts_layer_vary_number_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_layer_vary_number_f_shapes(1) = size(fortran_type_f%ts_layer_vary_number, 1)
  
  
  fortran_type_f%ts_n_totalcolumn_wfs = 0
  ts_n_totalcolumn_wfs_ptr => fortran_type_f%ts_n_totalcolumn_wfs
  transfer_struct_c%ts_n_totalcolumn_wfs = c_loc(ts_n_totalcolumn_wfs_ptr)
  inquire(iolength=transfer_struct_c%ts_n_totalcolumn_wfs_f_byte_size) fortran_type_f%ts_n_totalcolumn_wfs
#ifdef ifort
  transfer_struct_c%ts_n_totalcolumn_wfs_f_byte_size = transfer_struct_c%ts_n_totalcolumn_wfs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_n_surface_wfs = 0
  ts_n_surface_wfs_ptr => fortran_type_f%ts_n_surface_wfs
  transfer_struct_c%ts_n_surface_wfs = c_loc(ts_n_surface_wfs_ptr)
  inquire(iolength=transfer_struct_c%ts_n_surface_wfs_f_byte_size) fortran_type_f%ts_n_surface_wfs
#ifdef ifort
  transfer_struct_c%ts_n_surface_wfs_f_byte_size = transfer_struct_c%ts_n_surface_wfs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_n_sleave_wfs = 0
  ts_n_sleave_wfs_ptr => fortran_type_f%ts_n_sleave_wfs
  transfer_struct_c%ts_n_sleave_wfs = c_loc(ts_n_sleave_wfs_ptr)
  inquire(iolength=transfer_struct_c%ts_n_sleave_wfs_f_byte_size) fortran_type_f%ts_n_sleave_wfs
#ifdef ifort
  transfer_struct_c%ts_n_sleave_wfs_f_byte_size = transfer_struct_c%ts_n_sleave_wfs_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_columnwf_names = ''
  transfer_struct_c%ts_columnwf_names_f_len = len(fortran_type_f%ts_columnwf_names)
  transfer_struct_c%ts_columnwf_names_f_shapes(1) = size(fortran_type_f%ts_columnwf_names, 1)
  
  fortran_type_f%ts_profilewf_names = ''
  transfer_struct_c%ts_profilewf_names_f_len = len(fortran_type_f%ts_profilewf_names)
  transfer_struct_c%ts_profilewf_names_f_shapes(1) = size(fortran_type_f%ts_profilewf_names, 1)
  
  
end subroutine lidort_fixed_lincontrol_c_init_only

subroutine lidort_fixed_lincontrol_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_lincontrol_c_destroy

subroutine lidort_fixed_lincontrol_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f_from
  type(lidort_fixed_lincontrol), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_layer_vary_flag = fortran_type_f_from%ts_layer_vary_flag
  fortran_type_f_to%ts_layer_vary_number = fortran_type_f_from%ts_layer_vary_number
  fortran_type_f_to%ts_n_totalcolumn_wfs = fortran_type_f_from%ts_n_totalcolumn_wfs
  fortran_type_f_to%ts_n_surface_wfs = fortran_type_f_from%ts_n_surface_wfs
  fortran_type_f_to%ts_n_sleave_wfs = fortran_type_f_from%ts_n_sleave_wfs
  fortran_type_f_to%ts_columnwf_names = fortran_type_f_from%ts_columnwf_names
  fortran_type_f_to%ts_profilewf_names = fortran_type_f_from%ts_profilewf_names
  

end subroutine lidort_fixed_lincontrol_c_copy

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_linoptical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_linoptical

  type(lidort_fixed_linoptical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_linoptical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_linoptical_c_alloc_init

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_linoptical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_linoptical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  real(c_double), dimension(:,:), pointer :: ts_l_deltau_vert_input_ptr
  real(c_double), dimension(:,:), pointer :: ts_l_omega_total_input_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_l_phasmoms_total_input_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_l_phasfunc_input_up_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_l_phasfunc_input_dn_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_l_deltau_vert_input = 0_fpk
  ts_l_deltau_vert_input_ptr => fortran_type_f%ts_l_deltau_vert_input
  transfer_struct_c%ts_l_deltau_vert_input = c_loc(ts_l_deltau_vert_input_ptr(&
    lbound(fortran_type_f%ts_l_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,2)))
  inquire(iolength=transfer_struct_c%ts_l_deltau_vert_input_f_byte_size) fortran_type_f%ts_l_deltau_vert_input(&
    lbound(fortran_type_f%ts_l_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,2))
#ifdef ifort
  transfer_struct_c%ts_l_deltau_vert_input_f_byte_size = transfer_struct_c%ts_l_deltau_vert_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_deltau_vert_input_f_shapes(1) = size(fortran_type_f%ts_l_deltau_vert_input, 1)
  transfer_struct_c%ts_l_deltau_vert_input_f_shapes(2) = size(fortran_type_f%ts_l_deltau_vert_input, 2)
  
  
  fortran_type_f%ts_l_omega_total_input = 0_fpk
  ts_l_omega_total_input_ptr => fortran_type_f%ts_l_omega_total_input
  transfer_struct_c%ts_l_omega_total_input = c_loc(ts_l_omega_total_input_ptr(&
    lbound(fortran_type_f%ts_l_omega_total_input,1),&
    lbound(fortran_type_f%ts_l_omega_total_input,2)))
  inquire(iolength=transfer_struct_c%ts_l_omega_total_input_f_byte_size) fortran_type_f%ts_l_omega_total_input(&
    lbound(fortran_type_f%ts_l_omega_total_input,1),&
    lbound(fortran_type_f%ts_l_omega_total_input,2))
#ifdef ifort
  transfer_struct_c%ts_l_omega_total_input_f_byte_size = transfer_struct_c%ts_l_omega_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_omega_total_input_f_shapes(1) = size(fortran_type_f%ts_l_omega_total_input, 1)
  transfer_struct_c%ts_l_omega_total_input_f_shapes(2) = size(fortran_type_f%ts_l_omega_total_input, 2)
  
  
  fortran_type_f%ts_l_phasmoms_total_input = 0_fpk
  ts_l_phasmoms_total_input_ptr => fortran_type_f%ts_l_phasmoms_total_input
  transfer_struct_c%ts_l_phasmoms_total_input = c_loc(ts_l_phasmoms_total_input_ptr(&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,3)))
  inquire(iolength=transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size) fortran_type_f%ts_l_phasmoms_total_input(&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,3))
#ifdef ifort
  transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size = transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(1) = size(fortran_type_f%ts_l_phasmoms_total_input, 1)
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(2) = size(fortran_type_f%ts_l_phasmoms_total_input, 2)
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(3) = size(fortran_type_f%ts_l_phasmoms_total_input, 3)
  
  
  fortran_type_f%ts_l_phasfunc_input_up = 0_fpk
  ts_l_phasfunc_input_up_ptr => fortran_type_f%ts_l_phasfunc_input_up
  transfer_struct_c%ts_l_phasfunc_input_up = c_loc(ts_l_phasfunc_input_up_ptr(&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,1),&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,2),&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,3)))
  inquire(iolength=transfer_struct_c%ts_l_phasfunc_input_up_f_byte_size) fortran_type_f%ts_l_phasfunc_input_up(&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,1),&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,2),&
    lbound(fortran_type_f%ts_l_phasfunc_input_up,3))
#ifdef ifort
  transfer_struct_c%ts_l_phasfunc_input_up_f_byte_size = transfer_struct_c%ts_l_phasfunc_input_up_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_phasfunc_input_up_f_shapes(1) = size(fortran_type_f%ts_l_phasfunc_input_up, 1)
  transfer_struct_c%ts_l_phasfunc_input_up_f_shapes(2) = size(fortran_type_f%ts_l_phasfunc_input_up, 2)
  transfer_struct_c%ts_l_phasfunc_input_up_f_shapes(3) = size(fortran_type_f%ts_l_phasfunc_input_up, 3)
  
  
  fortran_type_f%ts_l_phasfunc_input_dn = 0_fpk
  ts_l_phasfunc_input_dn_ptr => fortran_type_f%ts_l_phasfunc_input_dn
  transfer_struct_c%ts_l_phasfunc_input_dn = c_loc(ts_l_phasfunc_input_dn_ptr(&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,1),&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,2),&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,3)))
  inquire(iolength=transfer_struct_c%ts_l_phasfunc_input_dn_f_byte_size) fortran_type_f%ts_l_phasfunc_input_dn(&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,1),&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,2),&
    lbound(fortran_type_f%ts_l_phasfunc_input_dn,3))
#ifdef ifort
  transfer_struct_c%ts_l_phasfunc_input_dn_f_byte_size = transfer_struct_c%ts_l_phasfunc_input_dn_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_phasfunc_input_dn_f_shapes(1) = size(fortran_type_f%ts_l_phasfunc_input_dn, 1)
  transfer_struct_c%ts_l_phasfunc_input_dn_f_shapes(2) = size(fortran_type_f%ts_l_phasfunc_input_dn, 2)
  transfer_struct_c%ts_l_phasfunc_input_dn_f_shapes(3) = size(fortran_type_f%ts_l_phasfunc_input_dn, 3)
  
  
end subroutine lidort_fixed_linoptical_c_init_only

subroutine lidort_fixed_linoptical_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_linoptical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_linoptical_c_destroy

subroutine lidort_fixed_linoptical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_linoptical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_linoptical), pointer :: fortran_type_f_from
  type(lidort_fixed_linoptical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_l_deltau_vert_input = fortran_type_f_from%ts_l_deltau_vert_input
  fortran_type_f_to%ts_l_omega_total_input = fortran_type_f_from%ts_l_omega_total_input
  fortran_type_f_to%ts_l_phasmoms_total_input = fortran_type_f_from%ts_l_phasmoms_total_input
  fortran_type_f_to%ts_l_phasfunc_input_up = fortran_type_f_from%ts_l_phasfunc_input_up
  fortran_type_f_to%ts_l_phasfunc_input_dn = fortran_type_f_from%ts_l_phasfunc_input_dn
  

end subroutine lidort_fixed_linoptical_c_copy

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lininputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lininputs

  type(lidort_fixed_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_lininputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_lininputs_c_alloc_init

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_lininputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  type(lidort_fixed_lincontrol), pointer :: cont_ptr
  type(lidort_fixed_linoptical), pointer :: optical_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  cont_ptr => fortran_type_f%cont
  transfer_struct_c%cont = c_loc(cont_ptr)
  inquire(iolength=transfer_struct_c%cont_f_byte_size) fortran_type_f%cont
#ifdef ifort
  transfer_struct_c%cont_f_byte_size = transfer_struct_c%cont_f_byte_size * 4
#endif
  
  
  optical_ptr => fortran_type_f%optical
  transfer_struct_c%optical = c_loc(optical_ptr)
  inquire(iolength=transfer_struct_c%optical_f_byte_size) fortran_type_f%optical
#ifdef ifort
  transfer_struct_c%optical_f_byte_size = transfer_struct_c%optical_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_lininputs_c_init_only

subroutine lidort_fixed_lininputs_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_lininputs_c_destroy

subroutine lidort_fixed_lininputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_lininputs), pointer :: fortran_type_f_from
  type(lidort_fixed_lininputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%cont = fortran_type_f_from%cont
  fortran_type_f_to%optical = fortran_type_f_from%optical
  

end subroutine lidort_fixed_lininputs_c_copy

! Links to type: "lidort_modified_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_lincontrol_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lincontrol

  type(lidort_modified_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lincontrol), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_lincontrol_c_alloc_init

! Links to type: "lidort_modified_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lincontrol), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_column_linearization_ptr
  logical(kind=4), pointer :: ts_do_profile_linearization_ptr
  logical(kind=4), pointer :: ts_do_atmos_linearization_ptr
  logical(kind=4), pointer :: ts_do_surface_linearization_ptr
  logical(kind=4), pointer :: ts_do_linearization_ptr
  logical(kind=4), pointer :: ts_do_simulation_only_ptr
  logical(kind=4), pointer :: ts_do_atmos_lbbf_ptr
  logical(kind=4), pointer :: ts_do_surface_lbbf_ptr
  logical(kind=4), pointer :: ts_do_sleave_wfs_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_do_column_linearization = .FALSE.
  ts_do_column_linearization_ptr => fortran_type_f%ts_do_column_linearization
  transfer_struct_c%ts_do_column_linearization = c_loc(ts_do_column_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_column_linearization_f_byte_size) fortran_type_f%ts_do_column_linearization
#ifdef ifort
  transfer_struct_c%ts_do_column_linearization_f_byte_size = transfer_struct_c%ts_do_column_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_profile_linearization = .FALSE.
  ts_do_profile_linearization_ptr => fortran_type_f%ts_do_profile_linearization
  transfer_struct_c%ts_do_profile_linearization = c_loc(ts_do_profile_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_profile_linearization_f_byte_size) fortran_type_f%ts_do_profile_linearization
#ifdef ifort
  transfer_struct_c%ts_do_profile_linearization_f_byte_size = transfer_struct_c%ts_do_profile_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_atmos_linearization = .FALSE.
  ts_do_atmos_linearization_ptr => fortran_type_f%ts_do_atmos_linearization
  transfer_struct_c%ts_do_atmos_linearization = c_loc(ts_do_atmos_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_atmos_linearization_f_byte_size) fortran_type_f%ts_do_atmos_linearization
#ifdef ifort
  transfer_struct_c%ts_do_atmos_linearization_f_byte_size = transfer_struct_c%ts_do_atmos_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_surface_linearization = .FALSE.
  ts_do_surface_linearization_ptr => fortran_type_f%ts_do_surface_linearization
  transfer_struct_c%ts_do_surface_linearization = c_loc(ts_do_surface_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_surface_linearization_f_byte_size) fortran_type_f%ts_do_surface_linearization
#ifdef ifort
  transfer_struct_c%ts_do_surface_linearization_f_byte_size = transfer_struct_c%ts_do_surface_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_linearization = .FALSE.
  ts_do_linearization_ptr => fortran_type_f%ts_do_linearization
  transfer_struct_c%ts_do_linearization = c_loc(ts_do_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_linearization_f_byte_size) fortran_type_f%ts_do_linearization
#ifdef ifort
  transfer_struct_c%ts_do_linearization_f_byte_size = transfer_struct_c%ts_do_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_simulation_only = .FALSE.
  ts_do_simulation_only_ptr => fortran_type_f%ts_do_simulation_only
  transfer_struct_c%ts_do_simulation_only = c_loc(ts_do_simulation_only_ptr)
  inquire(iolength=transfer_struct_c%ts_do_simulation_only_f_byte_size) fortran_type_f%ts_do_simulation_only
#ifdef ifort
  transfer_struct_c%ts_do_simulation_only_f_byte_size = transfer_struct_c%ts_do_simulation_only_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_atmos_lbbf = .FALSE.
  ts_do_atmos_lbbf_ptr => fortran_type_f%ts_do_atmos_lbbf
  transfer_struct_c%ts_do_atmos_lbbf = c_loc(ts_do_atmos_lbbf_ptr)
  inquire(iolength=transfer_struct_c%ts_do_atmos_lbbf_f_byte_size) fortran_type_f%ts_do_atmos_lbbf
#ifdef ifort
  transfer_struct_c%ts_do_atmos_lbbf_f_byte_size = transfer_struct_c%ts_do_atmos_lbbf_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_surface_lbbf = .FALSE.
  ts_do_surface_lbbf_ptr => fortran_type_f%ts_do_surface_lbbf
  transfer_struct_c%ts_do_surface_lbbf = c_loc(ts_do_surface_lbbf_ptr)
  inquire(iolength=transfer_struct_c%ts_do_surface_lbbf_f_byte_size) fortran_type_f%ts_do_surface_lbbf
#ifdef ifort
  transfer_struct_c%ts_do_surface_lbbf_f_byte_size = transfer_struct_c%ts_do_surface_lbbf_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sleave_wfs = .FALSE.
  ts_do_sleave_wfs_ptr => fortran_type_f%ts_do_sleave_wfs
  transfer_struct_c%ts_do_sleave_wfs = c_loc(ts_do_sleave_wfs_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sleave_wfs_f_byte_size) fortran_type_f%ts_do_sleave_wfs
#ifdef ifort
  transfer_struct_c%ts_do_sleave_wfs_f_byte_size = transfer_struct_c%ts_do_sleave_wfs_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_lincontrol_c_init_only

subroutine lidort_modified_lincontrol_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_lincontrol_c_destroy

subroutine lidort_modified_lincontrol_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_lincontrol), pointer :: fortran_type_f_from
  type(lidort_modified_lincontrol), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_column_linearization = fortran_type_f_from%ts_do_column_linearization
  fortran_type_f_to%ts_do_profile_linearization = fortran_type_f_from%ts_do_profile_linearization
  fortran_type_f_to%ts_do_atmos_linearization = fortran_type_f_from%ts_do_atmos_linearization
  fortran_type_f_to%ts_do_surface_linearization = fortran_type_f_from%ts_do_surface_linearization
  fortran_type_f_to%ts_do_linearization = fortran_type_f_from%ts_do_linearization
  fortran_type_f_to%ts_do_simulation_only = fortran_type_f_from%ts_do_simulation_only
  fortran_type_f_to%ts_do_atmos_lbbf = fortran_type_f_from%ts_do_atmos_lbbf
  fortran_type_f_to%ts_do_surface_lbbf = fortran_type_f_from%ts_do_surface_lbbf
  fortran_type_f_to%ts_do_sleave_wfs = fortran_type_f_from%ts_do_sleave_wfs
  

end subroutine lidort_modified_lincontrol_c_copy

! Links to type: "lidort_modified_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_lininputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lininputs

  type(lidort_modified_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_lininputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_lininputs_c_alloc_init

! Links to type: "lidort_modified_lininputs" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_lininputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  type(lidort_modified_lincontrol), pointer :: mcont_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  mcont_ptr => fortran_type_f%mcont
  transfer_struct_c%mcont = c_loc(mcont_ptr)
  inquire(iolength=transfer_struct_c%mcont_f_byte_size) fortran_type_f%mcont
#ifdef ifort
  transfer_struct_c%mcont_f_byte_size = transfer_struct_c%mcont_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_lininputs_c_init_only

subroutine lidort_modified_lininputs_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_lininputs_c_destroy

subroutine lidort_modified_lininputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_modified_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_lininputs), pointer :: fortran_type_f_from
  type(lidort_modified_lininputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%mcont = fortran_type_f_from%mcont
  

end subroutine lidort_modified_lininputs_c_copy

! Links to type: "lidort_linatmos" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linatmos_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linatmos

  type(lidort_linatmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linatmos_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linatmos_c_alloc_init

! Links to type: "lidort_linatmos" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linatmos_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m
  use lidort_pars_m

  type(lidort_linatmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_columnwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_meani_diffuse_colwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_flux_diffuse_colwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_dnmeani_direct_colwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_dnflux_direct_colwf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_profilewf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_meani_diffuse_profwf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_flux_diffuse_profwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_dnmeani_direct_profwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_dnflux_direct_profwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_abbwfs_jacobians_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_abbwfs_fluxes_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_albmed_user_profwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_trnmed_user_profwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_albmed_fluxes_profwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_trnmed_fluxes_profwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_transbeam_profwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_albmed_user_colwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_trnmed_user_colwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_albmed_fluxes_colwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_trnmed_fluxes_colwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_transbeam_colwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_planetary_transterm_profwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_planetary_sbterm_profwf_ptr
  real(c_double), dimension(:,:), pointer :: ts_planetary_transterm_colwf_ptr
  real(c_double), dimension(:), pointer :: ts_planetary_sbterm_colwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_lc_lostrans_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_lc_layer_mssts_ptr
  real(c_double), dimension(:,:), pointer :: ts_lc_surf_mssts_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_lp_lostrans_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_lp_layer_mssts_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_lp_surf_mssts_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_columnwf = 0_fpk
  ts_columnwf_ptr => fortran_type_f%ts_columnwf
  transfer_struct_c%ts_columnwf = c_loc(ts_columnwf_ptr(&
    lbound(fortran_type_f%ts_columnwf,1),&
    lbound(fortran_type_f%ts_columnwf,2),&
    lbound(fortran_type_f%ts_columnwf,3),&
    lbound(fortran_type_f%ts_columnwf,4)))
  inquire(iolength=transfer_struct_c%ts_columnwf_f_byte_size) fortran_type_f%ts_columnwf(&
    lbound(fortran_type_f%ts_columnwf,1),&
    lbound(fortran_type_f%ts_columnwf,2),&
    lbound(fortran_type_f%ts_columnwf,3),&
    lbound(fortran_type_f%ts_columnwf,4))
#ifdef ifort
  transfer_struct_c%ts_columnwf_f_byte_size = transfer_struct_c%ts_columnwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_columnwf_f_shapes(1) = size(fortran_type_f%ts_columnwf, 1)
  transfer_struct_c%ts_columnwf_f_shapes(2) = size(fortran_type_f%ts_columnwf, 2)
  transfer_struct_c%ts_columnwf_f_shapes(3) = size(fortran_type_f%ts_columnwf, 3)
  transfer_struct_c%ts_columnwf_f_shapes(4) = size(fortran_type_f%ts_columnwf, 4)
  
  
  fortran_type_f%ts_meani_diffuse_colwf = 0_fpk
  ts_meani_diffuse_colwf_ptr => fortran_type_f%ts_meani_diffuse_colwf
  transfer_struct_c%ts_meani_diffuse_colwf = c_loc(ts_meani_diffuse_colwf_ptr(&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,4)))
  inquire(iolength=transfer_struct_c%ts_meani_diffuse_colwf_f_byte_size) fortran_type_f%ts_meani_diffuse_colwf(&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_colwf,4))
#ifdef ifort
  transfer_struct_c%ts_meani_diffuse_colwf_f_byte_size = transfer_struct_c%ts_meani_diffuse_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_meani_diffuse_colwf_f_shapes(1) = size(fortran_type_f%ts_meani_diffuse_colwf, 1)
  transfer_struct_c%ts_meani_diffuse_colwf_f_shapes(2) = size(fortran_type_f%ts_meani_diffuse_colwf, 2)
  transfer_struct_c%ts_meani_diffuse_colwf_f_shapes(3) = size(fortran_type_f%ts_meani_diffuse_colwf, 3)
  transfer_struct_c%ts_meani_diffuse_colwf_f_shapes(4) = size(fortran_type_f%ts_meani_diffuse_colwf, 4)
  
  
  fortran_type_f%ts_flux_diffuse_colwf = 0_fpk
  ts_flux_diffuse_colwf_ptr => fortran_type_f%ts_flux_diffuse_colwf
  transfer_struct_c%ts_flux_diffuse_colwf = c_loc(ts_flux_diffuse_colwf_ptr(&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,4)))
  inquire(iolength=transfer_struct_c%ts_flux_diffuse_colwf_f_byte_size) fortran_type_f%ts_flux_diffuse_colwf(&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_colwf,4))
#ifdef ifort
  transfer_struct_c%ts_flux_diffuse_colwf_f_byte_size = transfer_struct_c%ts_flux_diffuse_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_diffuse_colwf_f_shapes(1) = size(fortran_type_f%ts_flux_diffuse_colwf, 1)
  transfer_struct_c%ts_flux_diffuse_colwf_f_shapes(2) = size(fortran_type_f%ts_flux_diffuse_colwf, 2)
  transfer_struct_c%ts_flux_diffuse_colwf_f_shapes(3) = size(fortran_type_f%ts_flux_diffuse_colwf, 3)
  transfer_struct_c%ts_flux_diffuse_colwf_f_shapes(4) = size(fortran_type_f%ts_flux_diffuse_colwf, 4)
  
  
  fortran_type_f%ts_dnmeani_direct_colwf = 0_fpk
  ts_dnmeani_direct_colwf_ptr => fortran_type_f%ts_dnmeani_direct_colwf
  transfer_struct_c%ts_dnmeani_direct_colwf = c_loc(ts_dnmeani_direct_colwf_ptr(&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,1),&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,2),&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,3)))
  inquire(iolength=transfer_struct_c%ts_dnmeani_direct_colwf_f_byte_size) fortran_type_f%ts_dnmeani_direct_colwf(&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,1),&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,2),&
    lbound(fortran_type_f%ts_dnmeani_direct_colwf,3))
#ifdef ifort
  transfer_struct_c%ts_dnmeani_direct_colwf_f_byte_size = transfer_struct_c%ts_dnmeani_direct_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnmeani_direct_colwf_f_shapes(1) = size(fortran_type_f%ts_dnmeani_direct_colwf, 1)
  transfer_struct_c%ts_dnmeani_direct_colwf_f_shapes(2) = size(fortran_type_f%ts_dnmeani_direct_colwf, 2)
  transfer_struct_c%ts_dnmeani_direct_colwf_f_shapes(3) = size(fortran_type_f%ts_dnmeani_direct_colwf, 3)
  
  
  fortran_type_f%ts_dnflux_direct_colwf = 0_fpk
  ts_dnflux_direct_colwf_ptr => fortran_type_f%ts_dnflux_direct_colwf
  transfer_struct_c%ts_dnflux_direct_colwf = c_loc(ts_dnflux_direct_colwf_ptr(&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,1),&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,2),&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,3)))
  inquire(iolength=transfer_struct_c%ts_dnflux_direct_colwf_f_byte_size) fortran_type_f%ts_dnflux_direct_colwf(&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,1),&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,2),&
    lbound(fortran_type_f%ts_dnflux_direct_colwf,3))
#ifdef ifort
  transfer_struct_c%ts_dnflux_direct_colwf_f_byte_size = transfer_struct_c%ts_dnflux_direct_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnflux_direct_colwf_f_shapes(1) = size(fortran_type_f%ts_dnflux_direct_colwf, 1)
  transfer_struct_c%ts_dnflux_direct_colwf_f_shapes(2) = size(fortran_type_f%ts_dnflux_direct_colwf, 2)
  transfer_struct_c%ts_dnflux_direct_colwf_f_shapes(3) = size(fortran_type_f%ts_dnflux_direct_colwf, 3)
  
  
  fortran_type_f%ts_profilewf = 0_fpk
  ts_profilewf_ptr => fortran_type_f%ts_profilewf
  transfer_struct_c%ts_profilewf = c_loc(ts_profilewf_ptr(&
    lbound(fortran_type_f%ts_profilewf,1),&
    lbound(fortran_type_f%ts_profilewf,2),&
    lbound(fortran_type_f%ts_profilewf,3),&
    lbound(fortran_type_f%ts_profilewf,4),&
    lbound(fortran_type_f%ts_profilewf,5)))
  inquire(iolength=transfer_struct_c%ts_profilewf_f_byte_size) fortran_type_f%ts_profilewf(&
    lbound(fortran_type_f%ts_profilewf,1),&
    lbound(fortran_type_f%ts_profilewf,2),&
    lbound(fortran_type_f%ts_profilewf,3),&
    lbound(fortran_type_f%ts_profilewf,4),&
    lbound(fortran_type_f%ts_profilewf,5))
#ifdef ifort
  transfer_struct_c%ts_profilewf_f_byte_size = transfer_struct_c%ts_profilewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_profilewf_f_shapes(1) = size(fortran_type_f%ts_profilewf, 1)
  transfer_struct_c%ts_profilewf_f_shapes(2) = size(fortran_type_f%ts_profilewf, 2)
  transfer_struct_c%ts_profilewf_f_shapes(3) = size(fortran_type_f%ts_profilewf, 3)
  transfer_struct_c%ts_profilewf_f_shapes(4) = size(fortran_type_f%ts_profilewf, 4)
  transfer_struct_c%ts_profilewf_f_shapes(5) = size(fortran_type_f%ts_profilewf, 5)
  
  
  fortran_type_f%ts_meani_diffuse_profwf = 0_fpk
  ts_meani_diffuse_profwf_ptr => fortran_type_f%ts_meani_diffuse_profwf
  transfer_struct_c%ts_meani_diffuse_profwf = c_loc(ts_meani_diffuse_profwf_ptr(&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,4),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,5)))
  inquire(iolength=transfer_struct_c%ts_meani_diffuse_profwf_f_byte_size) fortran_type_f%ts_meani_diffuse_profwf(&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,4),&
    lbound(fortran_type_f%ts_meani_diffuse_profwf,5))
#ifdef ifort
  transfer_struct_c%ts_meani_diffuse_profwf_f_byte_size = transfer_struct_c%ts_meani_diffuse_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_meani_diffuse_profwf_f_shapes(1) = size(fortran_type_f%ts_meani_diffuse_profwf, 1)
  transfer_struct_c%ts_meani_diffuse_profwf_f_shapes(2) = size(fortran_type_f%ts_meani_diffuse_profwf, 2)
  transfer_struct_c%ts_meani_diffuse_profwf_f_shapes(3) = size(fortran_type_f%ts_meani_diffuse_profwf, 3)
  transfer_struct_c%ts_meani_diffuse_profwf_f_shapes(4) = size(fortran_type_f%ts_meani_diffuse_profwf, 4)
  transfer_struct_c%ts_meani_diffuse_profwf_f_shapes(5) = size(fortran_type_f%ts_meani_diffuse_profwf, 5)
  
  
  fortran_type_f%ts_flux_diffuse_profwf = 0_fpk
  ts_flux_diffuse_profwf_ptr => fortran_type_f%ts_flux_diffuse_profwf
  transfer_struct_c%ts_flux_diffuse_profwf = c_loc(ts_flux_diffuse_profwf_ptr(&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,4),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,5)))
  inquire(iolength=transfer_struct_c%ts_flux_diffuse_profwf_f_byte_size) fortran_type_f%ts_flux_diffuse_profwf(&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,4),&
    lbound(fortran_type_f%ts_flux_diffuse_profwf,5))
#ifdef ifort
  transfer_struct_c%ts_flux_diffuse_profwf_f_byte_size = transfer_struct_c%ts_flux_diffuse_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_diffuse_profwf_f_shapes(1) = size(fortran_type_f%ts_flux_diffuse_profwf, 1)
  transfer_struct_c%ts_flux_diffuse_profwf_f_shapes(2) = size(fortran_type_f%ts_flux_diffuse_profwf, 2)
  transfer_struct_c%ts_flux_diffuse_profwf_f_shapes(3) = size(fortran_type_f%ts_flux_diffuse_profwf, 3)
  transfer_struct_c%ts_flux_diffuse_profwf_f_shapes(4) = size(fortran_type_f%ts_flux_diffuse_profwf, 4)
  transfer_struct_c%ts_flux_diffuse_profwf_f_shapes(5) = size(fortran_type_f%ts_flux_diffuse_profwf, 5)
  
  
  fortran_type_f%ts_dnmeani_direct_profwf = 0_fpk
  ts_dnmeani_direct_profwf_ptr => fortran_type_f%ts_dnmeani_direct_profwf
  transfer_struct_c%ts_dnmeani_direct_profwf = c_loc(ts_dnmeani_direct_profwf_ptr(&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,1),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,2),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,3),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,4)))
  inquire(iolength=transfer_struct_c%ts_dnmeani_direct_profwf_f_byte_size) fortran_type_f%ts_dnmeani_direct_profwf(&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,1),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,2),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,3),&
    lbound(fortran_type_f%ts_dnmeani_direct_profwf,4))
#ifdef ifort
  transfer_struct_c%ts_dnmeani_direct_profwf_f_byte_size = transfer_struct_c%ts_dnmeani_direct_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnmeani_direct_profwf_f_shapes(1) = size(fortran_type_f%ts_dnmeani_direct_profwf, 1)
  transfer_struct_c%ts_dnmeani_direct_profwf_f_shapes(2) = size(fortran_type_f%ts_dnmeani_direct_profwf, 2)
  transfer_struct_c%ts_dnmeani_direct_profwf_f_shapes(3) = size(fortran_type_f%ts_dnmeani_direct_profwf, 3)
  transfer_struct_c%ts_dnmeani_direct_profwf_f_shapes(4) = size(fortran_type_f%ts_dnmeani_direct_profwf, 4)
  
  
  fortran_type_f%ts_dnflux_direct_profwf = 0_fpk
  ts_dnflux_direct_profwf_ptr => fortran_type_f%ts_dnflux_direct_profwf
  transfer_struct_c%ts_dnflux_direct_profwf = c_loc(ts_dnflux_direct_profwf_ptr(&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,1),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,2),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,3),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,4)))
  inquire(iolength=transfer_struct_c%ts_dnflux_direct_profwf_f_byte_size) fortran_type_f%ts_dnflux_direct_profwf(&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,1),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,2),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,3),&
    lbound(fortran_type_f%ts_dnflux_direct_profwf,4))
#ifdef ifort
  transfer_struct_c%ts_dnflux_direct_profwf_f_byte_size = transfer_struct_c%ts_dnflux_direct_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnflux_direct_profwf_f_shapes(1) = size(fortran_type_f%ts_dnflux_direct_profwf, 1)
  transfer_struct_c%ts_dnflux_direct_profwf_f_shapes(2) = size(fortran_type_f%ts_dnflux_direct_profwf, 2)
  transfer_struct_c%ts_dnflux_direct_profwf_f_shapes(3) = size(fortran_type_f%ts_dnflux_direct_profwf, 3)
  transfer_struct_c%ts_dnflux_direct_profwf_f_shapes(4) = size(fortran_type_f%ts_dnflux_direct_profwf, 4)
  
  
  fortran_type_f%ts_abbwfs_jacobians = 0_fpk
  ts_abbwfs_jacobians_ptr => fortran_type_f%ts_abbwfs_jacobians
  transfer_struct_c%ts_abbwfs_jacobians = c_loc(ts_abbwfs_jacobians_ptr(&
    lbound(fortran_type_f%ts_abbwfs_jacobians,1),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,2),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,3),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,4)))
  inquire(iolength=transfer_struct_c%ts_abbwfs_jacobians_f_byte_size) fortran_type_f%ts_abbwfs_jacobians(&
    lbound(fortran_type_f%ts_abbwfs_jacobians,1),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,2),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,3),&
    lbound(fortran_type_f%ts_abbwfs_jacobians,4))
#ifdef ifort
  transfer_struct_c%ts_abbwfs_jacobians_f_byte_size = transfer_struct_c%ts_abbwfs_jacobians_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_abbwfs_jacobians_f_shapes(1) = size(fortran_type_f%ts_abbwfs_jacobians, 1)
  transfer_struct_c%ts_abbwfs_jacobians_f_shapes(2) = size(fortran_type_f%ts_abbwfs_jacobians, 2)
  transfer_struct_c%ts_abbwfs_jacobians_f_shapes(3) = size(fortran_type_f%ts_abbwfs_jacobians, 3)
  transfer_struct_c%ts_abbwfs_jacobians_f_shapes(4) = size(fortran_type_f%ts_abbwfs_jacobians, 4)
  
  
  fortran_type_f%ts_abbwfs_fluxes = 0_fpk
  ts_abbwfs_fluxes_ptr => fortran_type_f%ts_abbwfs_fluxes
  transfer_struct_c%ts_abbwfs_fluxes = c_loc(ts_abbwfs_fluxes_ptr(&
    lbound(fortran_type_f%ts_abbwfs_fluxes,1),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,2),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,3),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,4)))
  inquire(iolength=transfer_struct_c%ts_abbwfs_fluxes_f_byte_size) fortran_type_f%ts_abbwfs_fluxes(&
    lbound(fortran_type_f%ts_abbwfs_fluxes,1),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,2),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,3),&
    lbound(fortran_type_f%ts_abbwfs_fluxes,4))
#ifdef ifort
  transfer_struct_c%ts_abbwfs_fluxes_f_byte_size = transfer_struct_c%ts_abbwfs_fluxes_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_abbwfs_fluxes_f_shapes(1) = size(fortran_type_f%ts_abbwfs_fluxes, 1)
  transfer_struct_c%ts_abbwfs_fluxes_f_shapes(2) = size(fortran_type_f%ts_abbwfs_fluxes, 2)
  transfer_struct_c%ts_abbwfs_fluxes_f_shapes(3) = size(fortran_type_f%ts_abbwfs_fluxes, 3)
  transfer_struct_c%ts_abbwfs_fluxes_f_shapes(4) = size(fortran_type_f%ts_abbwfs_fluxes, 4)
  
  
  fortran_type_f%ts_albmed_user_profwf = 0_fpk
  ts_albmed_user_profwf_ptr => fortran_type_f%ts_albmed_user_profwf
  transfer_struct_c%ts_albmed_user_profwf = c_loc(ts_albmed_user_profwf_ptr(&
    lbound(fortran_type_f%ts_albmed_user_profwf,1),&
    lbound(fortran_type_f%ts_albmed_user_profwf,2),&
    lbound(fortran_type_f%ts_albmed_user_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_albmed_user_profwf_f_byte_size) fortran_type_f%ts_albmed_user_profwf(&
    lbound(fortran_type_f%ts_albmed_user_profwf,1),&
    lbound(fortran_type_f%ts_albmed_user_profwf,2),&
    lbound(fortran_type_f%ts_albmed_user_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_albmed_user_profwf_f_byte_size = transfer_struct_c%ts_albmed_user_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_user_profwf_f_shapes(1) = size(fortran_type_f%ts_albmed_user_profwf, 1)
  transfer_struct_c%ts_albmed_user_profwf_f_shapes(2) = size(fortran_type_f%ts_albmed_user_profwf, 2)
  transfer_struct_c%ts_albmed_user_profwf_f_shapes(3) = size(fortran_type_f%ts_albmed_user_profwf, 3)
  
  
  fortran_type_f%ts_trnmed_user_profwf = 0_fpk
  ts_trnmed_user_profwf_ptr => fortran_type_f%ts_trnmed_user_profwf
  transfer_struct_c%ts_trnmed_user_profwf = c_loc(ts_trnmed_user_profwf_ptr(&
    lbound(fortran_type_f%ts_trnmed_user_profwf,1),&
    lbound(fortran_type_f%ts_trnmed_user_profwf,2),&
    lbound(fortran_type_f%ts_trnmed_user_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_trnmed_user_profwf_f_byte_size) fortran_type_f%ts_trnmed_user_profwf(&
    lbound(fortran_type_f%ts_trnmed_user_profwf,1),&
    lbound(fortran_type_f%ts_trnmed_user_profwf,2),&
    lbound(fortran_type_f%ts_trnmed_user_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_trnmed_user_profwf_f_byte_size = transfer_struct_c%ts_trnmed_user_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_user_profwf_f_shapes(1) = size(fortran_type_f%ts_trnmed_user_profwf, 1)
  transfer_struct_c%ts_trnmed_user_profwf_f_shapes(2) = size(fortran_type_f%ts_trnmed_user_profwf, 2)
  transfer_struct_c%ts_trnmed_user_profwf_f_shapes(3) = size(fortran_type_f%ts_trnmed_user_profwf, 3)
  
  
  fortran_type_f%ts_albmed_fluxes_profwf = 0_fpk
  ts_albmed_fluxes_profwf_ptr => fortran_type_f%ts_albmed_fluxes_profwf
  transfer_struct_c%ts_albmed_fluxes_profwf = c_loc(ts_albmed_fluxes_profwf_ptr(&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,1),&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,2),&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_albmed_fluxes_profwf_f_byte_size) fortran_type_f%ts_albmed_fluxes_profwf(&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,1),&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,2),&
    lbound(fortran_type_f%ts_albmed_fluxes_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_albmed_fluxes_profwf_f_byte_size = transfer_struct_c%ts_albmed_fluxes_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_fluxes_profwf_f_shapes(1) = size(fortran_type_f%ts_albmed_fluxes_profwf, 1)
  transfer_struct_c%ts_albmed_fluxes_profwf_f_shapes(2) = size(fortran_type_f%ts_albmed_fluxes_profwf, 2)
  transfer_struct_c%ts_albmed_fluxes_profwf_f_shapes(3) = size(fortran_type_f%ts_albmed_fluxes_profwf, 3)
  
  
  fortran_type_f%ts_trnmed_fluxes_profwf = 0_fpk
  ts_trnmed_fluxes_profwf_ptr => fortran_type_f%ts_trnmed_fluxes_profwf
  transfer_struct_c%ts_trnmed_fluxes_profwf = c_loc(ts_trnmed_fluxes_profwf_ptr(&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,1),&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,2),&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_trnmed_fluxes_profwf_f_byte_size) fortran_type_f%ts_trnmed_fluxes_profwf(&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,1),&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,2),&
    lbound(fortran_type_f%ts_trnmed_fluxes_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_trnmed_fluxes_profwf_f_byte_size = transfer_struct_c%ts_trnmed_fluxes_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_fluxes_profwf_f_shapes(1) = size(fortran_type_f%ts_trnmed_fluxes_profwf, 1)
  transfer_struct_c%ts_trnmed_fluxes_profwf_f_shapes(2) = size(fortran_type_f%ts_trnmed_fluxes_profwf, 2)
  transfer_struct_c%ts_trnmed_fluxes_profwf_f_shapes(3) = size(fortran_type_f%ts_trnmed_fluxes_profwf, 3)
  
  
  fortran_type_f%ts_transbeam_profwf = 0_fpk
  ts_transbeam_profwf_ptr => fortran_type_f%ts_transbeam_profwf
  transfer_struct_c%ts_transbeam_profwf = c_loc(ts_transbeam_profwf_ptr(&
    lbound(fortran_type_f%ts_transbeam_profwf,1),&
    lbound(fortran_type_f%ts_transbeam_profwf,2),&
    lbound(fortran_type_f%ts_transbeam_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_transbeam_profwf_f_byte_size) fortran_type_f%ts_transbeam_profwf(&
    lbound(fortran_type_f%ts_transbeam_profwf,1),&
    lbound(fortran_type_f%ts_transbeam_profwf,2),&
    lbound(fortran_type_f%ts_transbeam_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_transbeam_profwf_f_byte_size = transfer_struct_c%ts_transbeam_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_transbeam_profwf_f_shapes(1) = size(fortran_type_f%ts_transbeam_profwf, 1)
  transfer_struct_c%ts_transbeam_profwf_f_shapes(2) = size(fortran_type_f%ts_transbeam_profwf, 2)
  transfer_struct_c%ts_transbeam_profwf_f_shapes(3) = size(fortran_type_f%ts_transbeam_profwf, 3)
  
  
  fortran_type_f%ts_albmed_user_colwf = 0_fpk
  ts_albmed_user_colwf_ptr => fortran_type_f%ts_albmed_user_colwf
  transfer_struct_c%ts_albmed_user_colwf = c_loc(ts_albmed_user_colwf_ptr(&
    lbound(fortran_type_f%ts_albmed_user_colwf,1),&
    lbound(fortran_type_f%ts_albmed_user_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_albmed_user_colwf_f_byte_size) fortran_type_f%ts_albmed_user_colwf(&
    lbound(fortran_type_f%ts_albmed_user_colwf,1),&
    lbound(fortran_type_f%ts_albmed_user_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_albmed_user_colwf_f_byte_size = transfer_struct_c%ts_albmed_user_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_user_colwf_f_shapes(1) = size(fortran_type_f%ts_albmed_user_colwf, 1)
  transfer_struct_c%ts_albmed_user_colwf_f_shapes(2) = size(fortran_type_f%ts_albmed_user_colwf, 2)
  
  
  fortran_type_f%ts_trnmed_user_colwf = 0_fpk
  ts_trnmed_user_colwf_ptr => fortran_type_f%ts_trnmed_user_colwf
  transfer_struct_c%ts_trnmed_user_colwf = c_loc(ts_trnmed_user_colwf_ptr(&
    lbound(fortran_type_f%ts_trnmed_user_colwf,1),&
    lbound(fortran_type_f%ts_trnmed_user_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_trnmed_user_colwf_f_byte_size) fortran_type_f%ts_trnmed_user_colwf(&
    lbound(fortran_type_f%ts_trnmed_user_colwf,1),&
    lbound(fortran_type_f%ts_trnmed_user_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_trnmed_user_colwf_f_byte_size = transfer_struct_c%ts_trnmed_user_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_user_colwf_f_shapes(1) = size(fortran_type_f%ts_trnmed_user_colwf, 1)
  transfer_struct_c%ts_trnmed_user_colwf_f_shapes(2) = size(fortran_type_f%ts_trnmed_user_colwf, 2)
  
  
  fortran_type_f%ts_albmed_fluxes_colwf = 0_fpk
  ts_albmed_fluxes_colwf_ptr => fortran_type_f%ts_albmed_fluxes_colwf
  transfer_struct_c%ts_albmed_fluxes_colwf = c_loc(ts_albmed_fluxes_colwf_ptr(&
    lbound(fortran_type_f%ts_albmed_fluxes_colwf,1),&
    lbound(fortran_type_f%ts_albmed_fluxes_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_albmed_fluxes_colwf_f_byte_size) fortran_type_f%ts_albmed_fluxes_colwf(&
    lbound(fortran_type_f%ts_albmed_fluxes_colwf,1),&
    lbound(fortran_type_f%ts_albmed_fluxes_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_albmed_fluxes_colwf_f_byte_size = transfer_struct_c%ts_albmed_fluxes_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_fluxes_colwf_f_shapes(1) = size(fortran_type_f%ts_albmed_fluxes_colwf, 1)
  transfer_struct_c%ts_albmed_fluxes_colwf_f_shapes(2) = size(fortran_type_f%ts_albmed_fluxes_colwf, 2)
  
  
  fortran_type_f%ts_trnmed_fluxes_colwf = 0_fpk
  ts_trnmed_fluxes_colwf_ptr => fortran_type_f%ts_trnmed_fluxes_colwf
  transfer_struct_c%ts_trnmed_fluxes_colwf = c_loc(ts_trnmed_fluxes_colwf_ptr(&
    lbound(fortran_type_f%ts_trnmed_fluxes_colwf,1),&
    lbound(fortran_type_f%ts_trnmed_fluxes_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_trnmed_fluxes_colwf_f_byte_size) fortran_type_f%ts_trnmed_fluxes_colwf(&
    lbound(fortran_type_f%ts_trnmed_fluxes_colwf,1),&
    lbound(fortran_type_f%ts_trnmed_fluxes_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_trnmed_fluxes_colwf_f_byte_size = transfer_struct_c%ts_trnmed_fluxes_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_fluxes_colwf_f_shapes(1) = size(fortran_type_f%ts_trnmed_fluxes_colwf, 1)
  transfer_struct_c%ts_trnmed_fluxes_colwf_f_shapes(2) = size(fortran_type_f%ts_trnmed_fluxes_colwf, 2)
  
  
  fortran_type_f%ts_transbeam_colwf = 0_fpk
  ts_transbeam_colwf_ptr => fortran_type_f%ts_transbeam_colwf
  transfer_struct_c%ts_transbeam_colwf = c_loc(ts_transbeam_colwf_ptr(&
    lbound(fortran_type_f%ts_transbeam_colwf,1),&
    lbound(fortran_type_f%ts_transbeam_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_transbeam_colwf_f_byte_size) fortran_type_f%ts_transbeam_colwf(&
    lbound(fortran_type_f%ts_transbeam_colwf,1),&
    lbound(fortran_type_f%ts_transbeam_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_transbeam_colwf_f_byte_size = transfer_struct_c%ts_transbeam_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_transbeam_colwf_f_shapes(1) = size(fortran_type_f%ts_transbeam_colwf, 1)
  transfer_struct_c%ts_transbeam_colwf_f_shapes(2) = size(fortran_type_f%ts_transbeam_colwf, 2)
  
  
  fortran_type_f%ts_planetary_transterm_profwf = 0_fpk
  ts_planetary_transterm_profwf_ptr => fortran_type_f%ts_planetary_transterm_profwf
  transfer_struct_c%ts_planetary_transterm_profwf = c_loc(ts_planetary_transterm_profwf_ptr(&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,1),&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,2),&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,3)))
  inquire(iolength=transfer_struct_c%ts_planetary_transterm_profwf_f_byte_size) fortran_type_f%ts_planetary_transterm_profwf(&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,1),&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,2),&
    lbound(fortran_type_f%ts_planetary_transterm_profwf,3))
#ifdef ifort
  transfer_struct_c%ts_planetary_transterm_profwf_f_byte_size = transfer_struct_c%ts_planetary_transterm_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_planetary_transterm_profwf_f_shapes(1) = size(fortran_type_f%ts_planetary_transterm_profwf, 1)
  transfer_struct_c%ts_planetary_transterm_profwf_f_shapes(2) = size(fortran_type_f%ts_planetary_transterm_profwf, 2)
  transfer_struct_c%ts_planetary_transterm_profwf_f_shapes(3) = size(fortran_type_f%ts_planetary_transterm_profwf, 3)
  
  
  fortran_type_f%ts_planetary_sbterm_profwf = 0_fpk
  ts_planetary_sbterm_profwf_ptr => fortran_type_f%ts_planetary_sbterm_profwf
  transfer_struct_c%ts_planetary_sbterm_profwf = c_loc(ts_planetary_sbterm_profwf_ptr(&
    lbound(fortran_type_f%ts_planetary_sbterm_profwf,1),&
    lbound(fortran_type_f%ts_planetary_sbterm_profwf,2)))
  inquire(iolength=transfer_struct_c%ts_planetary_sbterm_profwf_f_byte_size) fortran_type_f%ts_planetary_sbterm_profwf(&
    lbound(fortran_type_f%ts_planetary_sbterm_profwf,1),&
    lbound(fortran_type_f%ts_planetary_sbterm_profwf,2))
#ifdef ifort
  transfer_struct_c%ts_planetary_sbterm_profwf_f_byte_size = transfer_struct_c%ts_planetary_sbterm_profwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_planetary_sbterm_profwf_f_shapes(1) = size(fortran_type_f%ts_planetary_sbterm_profwf, 1)
  transfer_struct_c%ts_planetary_sbterm_profwf_f_shapes(2) = size(fortran_type_f%ts_planetary_sbterm_profwf, 2)
  
  
  fortran_type_f%ts_planetary_transterm_colwf = 0_fpk
  ts_planetary_transterm_colwf_ptr => fortran_type_f%ts_planetary_transterm_colwf
  transfer_struct_c%ts_planetary_transterm_colwf = c_loc(ts_planetary_transterm_colwf_ptr(&
    lbound(fortran_type_f%ts_planetary_transterm_colwf,1),&
    lbound(fortran_type_f%ts_planetary_transterm_colwf,2)))
  inquire(iolength=transfer_struct_c%ts_planetary_transterm_colwf_f_byte_size) fortran_type_f%ts_planetary_transterm_colwf(&
    lbound(fortran_type_f%ts_planetary_transterm_colwf,1),&
    lbound(fortran_type_f%ts_planetary_transterm_colwf,2))
#ifdef ifort
  transfer_struct_c%ts_planetary_transterm_colwf_f_byte_size = transfer_struct_c%ts_planetary_transterm_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_planetary_transterm_colwf_f_shapes(1) = size(fortran_type_f%ts_planetary_transterm_colwf, 1)
  transfer_struct_c%ts_planetary_transterm_colwf_f_shapes(2) = size(fortran_type_f%ts_planetary_transterm_colwf, 2)
  
  
  fortran_type_f%ts_planetary_sbterm_colwf = 0_fpk
  ts_planetary_sbterm_colwf_ptr => fortran_type_f%ts_planetary_sbterm_colwf
  transfer_struct_c%ts_planetary_sbterm_colwf = c_loc(ts_planetary_sbterm_colwf_ptr(&
    lbound(fortran_type_f%ts_planetary_sbterm_colwf,1)))
  inquire(iolength=transfer_struct_c%ts_planetary_sbterm_colwf_f_byte_size) fortran_type_f%ts_planetary_sbterm_colwf(&
    lbound(fortran_type_f%ts_planetary_sbterm_colwf,1))
#ifdef ifort
  transfer_struct_c%ts_planetary_sbterm_colwf_f_byte_size = transfer_struct_c%ts_planetary_sbterm_colwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_planetary_sbterm_colwf_f_shapes(1) = size(fortran_type_f%ts_planetary_sbterm_colwf, 1)
  
  
  fortran_type_f%ts_lc_lostrans = 0_fpk
  ts_lc_lostrans_ptr => fortran_type_f%ts_lc_lostrans
  transfer_struct_c%ts_lc_lostrans = c_loc(ts_lc_lostrans_ptr(&
    lbound(fortran_type_f%ts_lc_lostrans,1),&
    lbound(fortran_type_f%ts_lc_lostrans,2),&
    lbound(fortran_type_f%ts_lc_lostrans,3)))
  inquire(iolength=transfer_struct_c%ts_lc_lostrans_f_byte_size) fortran_type_f%ts_lc_lostrans(&
    lbound(fortran_type_f%ts_lc_lostrans,1),&
    lbound(fortran_type_f%ts_lc_lostrans,2),&
    lbound(fortran_type_f%ts_lc_lostrans,3))
#ifdef ifort
  transfer_struct_c%ts_lc_lostrans_f_byte_size = transfer_struct_c%ts_lc_lostrans_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lc_lostrans_f_shapes(1) = size(fortran_type_f%ts_lc_lostrans, 1)
  transfer_struct_c%ts_lc_lostrans_f_shapes(2) = size(fortran_type_f%ts_lc_lostrans, 2)
  transfer_struct_c%ts_lc_lostrans_f_shapes(3) = size(fortran_type_f%ts_lc_lostrans, 3)
  
  
  fortran_type_f%ts_lc_layer_mssts = 0_fpk
  ts_lc_layer_mssts_ptr => fortran_type_f%ts_lc_layer_mssts
  transfer_struct_c%ts_lc_layer_mssts = c_loc(ts_lc_layer_mssts_ptr(&
    lbound(fortran_type_f%ts_lc_layer_mssts,1),&
    lbound(fortran_type_f%ts_lc_layer_mssts,2),&
    lbound(fortran_type_f%ts_lc_layer_mssts,3)))
  inquire(iolength=transfer_struct_c%ts_lc_layer_mssts_f_byte_size) fortran_type_f%ts_lc_layer_mssts(&
    lbound(fortran_type_f%ts_lc_layer_mssts,1),&
    lbound(fortran_type_f%ts_lc_layer_mssts,2),&
    lbound(fortran_type_f%ts_lc_layer_mssts,3))
#ifdef ifort
  transfer_struct_c%ts_lc_layer_mssts_f_byte_size = transfer_struct_c%ts_lc_layer_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lc_layer_mssts_f_shapes(1) = size(fortran_type_f%ts_lc_layer_mssts, 1)
  transfer_struct_c%ts_lc_layer_mssts_f_shapes(2) = size(fortran_type_f%ts_lc_layer_mssts, 2)
  transfer_struct_c%ts_lc_layer_mssts_f_shapes(3) = size(fortran_type_f%ts_lc_layer_mssts, 3)
  
  
  fortran_type_f%ts_lc_surf_mssts = 0_fpk
  ts_lc_surf_mssts_ptr => fortran_type_f%ts_lc_surf_mssts
  transfer_struct_c%ts_lc_surf_mssts = c_loc(ts_lc_surf_mssts_ptr(&
    lbound(fortran_type_f%ts_lc_surf_mssts,1),&
    lbound(fortran_type_f%ts_lc_surf_mssts,2)))
  inquire(iolength=transfer_struct_c%ts_lc_surf_mssts_f_byte_size) fortran_type_f%ts_lc_surf_mssts(&
    lbound(fortran_type_f%ts_lc_surf_mssts,1),&
    lbound(fortran_type_f%ts_lc_surf_mssts,2))
#ifdef ifort
  transfer_struct_c%ts_lc_surf_mssts_f_byte_size = transfer_struct_c%ts_lc_surf_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lc_surf_mssts_f_shapes(1) = size(fortran_type_f%ts_lc_surf_mssts, 1)
  transfer_struct_c%ts_lc_surf_mssts_f_shapes(2) = size(fortran_type_f%ts_lc_surf_mssts, 2)
  
  
  fortran_type_f%ts_lp_lostrans = 0_fpk
  ts_lp_lostrans_ptr => fortran_type_f%ts_lp_lostrans
  transfer_struct_c%ts_lp_lostrans = c_loc(ts_lp_lostrans_ptr(&
    lbound(fortran_type_f%ts_lp_lostrans,1),&
    lbound(fortran_type_f%ts_lp_lostrans,2),&
    lbound(fortran_type_f%ts_lp_lostrans,3)))
  inquire(iolength=transfer_struct_c%ts_lp_lostrans_f_byte_size) fortran_type_f%ts_lp_lostrans(&
    lbound(fortran_type_f%ts_lp_lostrans,1),&
    lbound(fortran_type_f%ts_lp_lostrans,2),&
    lbound(fortran_type_f%ts_lp_lostrans,3))
#ifdef ifort
  transfer_struct_c%ts_lp_lostrans_f_byte_size = transfer_struct_c%ts_lp_lostrans_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lp_lostrans_f_shapes(1) = size(fortran_type_f%ts_lp_lostrans, 1)
  transfer_struct_c%ts_lp_lostrans_f_shapes(2) = size(fortran_type_f%ts_lp_lostrans, 2)
  transfer_struct_c%ts_lp_lostrans_f_shapes(3) = size(fortran_type_f%ts_lp_lostrans, 3)
  
  
  fortran_type_f%ts_lp_layer_mssts = 0_fpk
  ts_lp_layer_mssts_ptr => fortran_type_f%ts_lp_layer_mssts
  transfer_struct_c%ts_lp_layer_mssts = c_loc(ts_lp_layer_mssts_ptr(&
    lbound(fortran_type_f%ts_lp_layer_mssts,1),&
    lbound(fortran_type_f%ts_lp_layer_mssts,2),&
    lbound(fortran_type_f%ts_lp_layer_mssts,3),&
    lbound(fortran_type_f%ts_lp_layer_mssts,4)))
  inquire(iolength=transfer_struct_c%ts_lp_layer_mssts_f_byte_size) fortran_type_f%ts_lp_layer_mssts(&
    lbound(fortran_type_f%ts_lp_layer_mssts,1),&
    lbound(fortran_type_f%ts_lp_layer_mssts,2),&
    lbound(fortran_type_f%ts_lp_layer_mssts,3),&
    lbound(fortran_type_f%ts_lp_layer_mssts,4))
#ifdef ifort
  transfer_struct_c%ts_lp_layer_mssts_f_byte_size = transfer_struct_c%ts_lp_layer_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lp_layer_mssts_f_shapes(1) = size(fortran_type_f%ts_lp_layer_mssts, 1)
  transfer_struct_c%ts_lp_layer_mssts_f_shapes(2) = size(fortran_type_f%ts_lp_layer_mssts, 2)
  transfer_struct_c%ts_lp_layer_mssts_f_shapes(3) = size(fortran_type_f%ts_lp_layer_mssts, 3)
  transfer_struct_c%ts_lp_layer_mssts_f_shapes(4) = size(fortran_type_f%ts_lp_layer_mssts, 4)
  
  
  fortran_type_f%ts_lp_surf_mssts = 0_fpk
  ts_lp_surf_mssts_ptr => fortran_type_f%ts_lp_surf_mssts
  transfer_struct_c%ts_lp_surf_mssts = c_loc(ts_lp_surf_mssts_ptr(&
    lbound(fortran_type_f%ts_lp_surf_mssts,1),&
    lbound(fortran_type_f%ts_lp_surf_mssts,2),&
    lbound(fortran_type_f%ts_lp_surf_mssts,3)))
  inquire(iolength=transfer_struct_c%ts_lp_surf_mssts_f_byte_size) fortran_type_f%ts_lp_surf_mssts(&
    lbound(fortran_type_f%ts_lp_surf_mssts,1),&
    lbound(fortran_type_f%ts_lp_surf_mssts,2),&
    lbound(fortran_type_f%ts_lp_surf_mssts,3))
#ifdef ifort
  transfer_struct_c%ts_lp_surf_mssts_f_byte_size = transfer_struct_c%ts_lp_surf_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lp_surf_mssts_f_shapes(1) = size(fortran_type_f%ts_lp_surf_mssts, 1)
  transfer_struct_c%ts_lp_surf_mssts_f_shapes(2) = size(fortran_type_f%ts_lp_surf_mssts, 2)
  transfer_struct_c%ts_lp_surf_mssts_f_shapes(3) = size(fortran_type_f%ts_lp_surf_mssts, 3)
  
  
end subroutine lidort_linatmos_c_init_only

subroutine lidort_linatmos_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linatmos

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linatmos_c_destroy

subroutine lidort_linatmos_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linatmos

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linatmos), pointer :: fortran_type_f_from
  type(lidort_linatmos), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_columnwf = fortran_type_f_from%ts_columnwf
  fortran_type_f_to%ts_meani_diffuse_colwf = fortran_type_f_from%ts_meani_diffuse_colwf
  fortran_type_f_to%ts_flux_diffuse_colwf = fortran_type_f_from%ts_flux_diffuse_colwf
  fortran_type_f_to%ts_dnmeani_direct_colwf = fortran_type_f_from%ts_dnmeani_direct_colwf
  fortran_type_f_to%ts_dnflux_direct_colwf = fortran_type_f_from%ts_dnflux_direct_colwf
  fortran_type_f_to%ts_profilewf = fortran_type_f_from%ts_profilewf
  fortran_type_f_to%ts_meani_diffuse_profwf = fortran_type_f_from%ts_meani_diffuse_profwf
  fortran_type_f_to%ts_flux_diffuse_profwf = fortran_type_f_from%ts_flux_diffuse_profwf
  fortran_type_f_to%ts_dnmeani_direct_profwf = fortran_type_f_from%ts_dnmeani_direct_profwf
  fortran_type_f_to%ts_dnflux_direct_profwf = fortran_type_f_from%ts_dnflux_direct_profwf
  fortran_type_f_to%ts_abbwfs_jacobians = fortran_type_f_from%ts_abbwfs_jacobians
  fortran_type_f_to%ts_abbwfs_fluxes = fortran_type_f_from%ts_abbwfs_fluxes
  fortran_type_f_to%ts_albmed_user_profwf = fortran_type_f_from%ts_albmed_user_profwf
  fortran_type_f_to%ts_trnmed_user_profwf = fortran_type_f_from%ts_trnmed_user_profwf
  fortran_type_f_to%ts_albmed_fluxes_profwf = fortran_type_f_from%ts_albmed_fluxes_profwf
  fortran_type_f_to%ts_trnmed_fluxes_profwf = fortran_type_f_from%ts_trnmed_fluxes_profwf
  fortran_type_f_to%ts_transbeam_profwf = fortran_type_f_from%ts_transbeam_profwf
  fortran_type_f_to%ts_albmed_user_colwf = fortran_type_f_from%ts_albmed_user_colwf
  fortran_type_f_to%ts_trnmed_user_colwf = fortran_type_f_from%ts_trnmed_user_colwf
  fortran_type_f_to%ts_albmed_fluxes_colwf = fortran_type_f_from%ts_albmed_fluxes_colwf
  fortran_type_f_to%ts_trnmed_fluxes_colwf = fortran_type_f_from%ts_trnmed_fluxes_colwf
  fortran_type_f_to%ts_transbeam_colwf = fortran_type_f_from%ts_transbeam_colwf
  fortran_type_f_to%ts_planetary_transterm_profwf = fortran_type_f_from%ts_planetary_transterm_profwf
  fortran_type_f_to%ts_planetary_sbterm_profwf = fortran_type_f_from%ts_planetary_sbterm_profwf
  fortran_type_f_to%ts_planetary_transterm_colwf = fortran_type_f_from%ts_planetary_transterm_colwf
  fortran_type_f_to%ts_planetary_sbterm_colwf = fortran_type_f_from%ts_planetary_sbterm_colwf
  fortran_type_f_to%ts_lc_lostrans = fortran_type_f_from%ts_lc_lostrans
  fortran_type_f_to%ts_lc_layer_mssts = fortran_type_f_from%ts_lc_layer_mssts
  fortran_type_f_to%ts_lc_surf_mssts = fortran_type_f_from%ts_lc_surf_mssts
  fortran_type_f_to%ts_lp_lostrans = fortran_type_f_from%ts_lp_lostrans
  fortran_type_f_to%ts_lp_layer_mssts = fortran_type_f_from%ts_lp_layer_mssts
  fortran_type_f_to%ts_lp_surf_mssts = fortran_type_f_from%ts_lp_surf_mssts
  

end subroutine lidort_linatmos_c_copy

! Links to type: "lidort_linsurf" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linsurf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linsurf

  type(lidort_linsurf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsurf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsurf_c_alloc_init

! Links to type: "lidort_linsurf" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linsurf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m
  use lidort_pars_m

  type(lidort_linsurf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_surfacewf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_meani_diffuse_surfwf_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_flux_diffuse_surfwf_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_sbbwfs_jacobians_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_sbbwfs_fluxes_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_ls_layer_mssts_ptr
  real(c_double), dimension(:,:), pointer :: ts_ls_surf_mssts_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_surfacewf = 0_fpk
  ts_surfacewf_ptr => fortran_type_f%ts_surfacewf
  transfer_struct_c%ts_surfacewf = c_loc(ts_surfacewf_ptr(&
    lbound(fortran_type_f%ts_surfacewf,1),&
    lbound(fortran_type_f%ts_surfacewf,2),&
    lbound(fortran_type_f%ts_surfacewf,3),&
    lbound(fortran_type_f%ts_surfacewf,4)))
  inquire(iolength=transfer_struct_c%ts_surfacewf_f_byte_size) fortran_type_f%ts_surfacewf(&
    lbound(fortran_type_f%ts_surfacewf,1),&
    lbound(fortran_type_f%ts_surfacewf,2),&
    lbound(fortran_type_f%ts_surfacewf,3),&
    lbound(fortran_type_f%ts_surfacewf,4))
#ifdef ifort
  transfer_struct_c%ts_surfacewf_f_byte_size = transfer_struct_c%ts_surfacewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_surfacewf_f_shapes(1) = size(fortran_type_f%ts_surfacewf, 1)
  transfer_struct_c%ts_surfacewf_f_shapes(2) = size(fortran_type_f%ts_surfacewf, 2)
  transfer_struct_c%ts_surfacewf_f_shapes(3) = size(fortran_type_f%ts_surfacewf, 3)
  transfer_struct_c%ts_surfacewf_f_shapes(4) = size(fortran_type_f%ts_surfacewf, 4)
  
  
  fortran_type_f%ts_meani_diffuse_surfwf = 0_fpk
  ts_meani_diffuse_surfwf_ptr => fortran_type_f%ts_meani_diffuse_surfwf
  transfer_struct_c%ts_meani_diffuse_surfwf = c_loc(ts_meani_diffuse_surfwf_ptr(&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,4)))
  inquire(iolength=transfer_struct_c%ts_meani_diffuse_surfwf_f_byte_size) fortran_type_f%ts_meani_diffuse_surfwf(&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,1),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,2),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,3),&
    lbound(fortran_type_f%ts_meani_diffuse_surfwf,4))
#ifdef ifort
  transfer_struct_c%ts_meani_diffuse_surfwf_f_byte_size = transfer_struct_c%ts_meani_diffuse_surfwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_meani_diffuse_surfwf_f_shapes(1) = size(fortran_type_f%ts_meani_diffuse_surfwf, 1)
  transfer_struct_c%ts_meani_diffuse_surfwf_f_shapes(2) = size(fortran_type_f%ts_meani_diffuse_surfwf, 2)
  transfer_struct_c%ts_meani_diffuse_surfwf_f_shapes(3) = size(fortran_type_f%ts_meani_diffuse_surfwf, 3)
  transfer_struct_c%ts_meani_diffuse_surfwf_f_shapes(4) = size(fortran_type_f%ts_meani_diffuse_surfwf, 4)
  
  
  fortran_type_f%ts_flux_diffuse_surfwf = 0_fpk
  ts_flux_diffuse_surfwf_ptr => fortran_type_f%ts_flux_diffuse_surfwf
  transfer_struct_c%ts_flux_diffuse_surfwf = c_loc(ts_flux_diffuse_surfwf_ptr(&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,4)))
  inquire(iolength=transfer_struct_c%ts_flux_diffuse_surfwf_f_byte_size) fortran_type_f%ts_flux_diffuse_surfwf(&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,1),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,2),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,3),&
    lbound(fortran_type_f%ts_flux_diffuse_surfwf,4))
#ifdef ifort
  transfer_struct_c%ts_flux_diffuse_surfwf_f_byte_size = transfer_struct_c%ts_flux_diffuse_surfwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_diffuse_surfwf_f_shapes(1) = size(fortran_type_f%ts_flux_diffuse_surfwf, 1)
  transfer_struct_c%ts_flux_diffuse_surfwf_f_shapes(2) = size(fortran_type_f%ts_flux_diffuse_surfwf, 2)
  transfer_struct_c%ts_flux_diffuse_surfwf_f_shapes(3) = size(fortran_type_f%ts_flux_diffuse_surfwf, 3)
  transfer_struct_c%ts_flux_diffuse_surfwf_f_shapes(4) = size(fortran_type_f%ts_flux_diffuse_surfwf, 4)
  
  
  fortran_type_f%ts_sbbwfs_jacobians = 0_fpk
  ts_sbbwfs_jacobians_ptr => fortran_type_f%ts_sbbwfs_jacobians
  transfer_struct_c%ts_sbbwfs_jacobians = c_loc(ts_sbbwfs_jacobians_ptr(&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,1),&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,2),&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,3)))
  inquire(iolength=transfer_struct_c%ts_sbbwfs_jacobians_f_byte_size) fortran_type_f%ts_sbbwfs_jacobians(&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,1),&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,2),&
    lbound(fortran_type_f%ts_sbbwfs_jacobians,3))
#ifdef ifort
  transfer_struct_c%ts_sbbwfs_jacobians_f_byte_size = transfer_struct_c%ts_sbbwfs_jacobians_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_sbbwfs_jacobians_f_shapes(1) = size(fortran_type_f%ts_sbbwfs_jacobians, 1)
  transfer_struct_c%ts_sbbwfs_jacobians_f_shapes(2) = size(fortran_type_f%ts_sbbwfs_jacobians, 2)
  transfer_struct_c%ts_sbbwfs_jacobians_f_shapes(3) = size(fortran_type_f%ts_sbbwfs_jacobians, 3)
  
  
  fortran_type_f%ts_sbbwfs_fluxes = 0_fpk
  ts_sbbwfs_fluxes_ptr => fortran_type_f%ts_sbbwfs_fluxes
  transfer_struct_c%ts_sbbwfs_fluxes = c_loc(ts_sbbwfs_fluxes_ptr(&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,1),&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,2),&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,3)))
  inquire(iolength=transfer_struct_c%ts_sbbwfs_fluxes_f_byte_size) fortran_type_f%ts_sbbwfs_fluxes(&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,1),&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,2),&
    lbound(fortran_type_f%ts_sbbwfs_fluxes,3))
#ifdef ifort
  transfer_struct_c%ts_sbbwfs_fluxes_f_byte_size = transfer_struct_c%ts_sbbwfs_fluxes_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_sbbwfs_fluxes_f_shapes(1) = size(fortran_type_f%ts_sbbwfs_fluxes, 1)
  transfer_struct_c%ts_sbbwfs_fluxes_f_shapes(2) = size(fortran_type_f%ts_sbbwfs_fluxes, 2)
  transfer_struct_c%ts_sbbwfs_fluxes_f_shapes(3) = size(fortran_type_f%ts_sbbwfs_fluxes, 3)
  
  
  fortran_type_f%ts_ls_layer_mssts = 0_fpk
  ts_ls_layer_mssts_ptr => fortran_type_f%ts_ls_layer_mssts
  transfer_struct_c%ts_ls_layer_mssts = c_loc(ts_ls_layer_mssts_ptr(&
    lbound(fortran_type_f%ts_ls_layer_mssts,1),&
    lbound(fortran_type_f%ts_ls_layer_mssts,2),&
    lbound(fortran_type_f%ts_ls_layer_mssts,3)))
  inquire(iolength=transfer_struct_c%ts_ls_layer_mssts_f_byte_size) fortran_type_f%ts_ls_layer_mssts(&
    lbound(fortran_type_f%ts_ls_layer_mssts,1),&
    lbound(fortran_type_f%ts_ls_layer_mssts,2),&
    lbound(fortran_type_f%ts_ls_layer_mssts,3))
#ifdef ifort
  transfer_struct_c%ts_ls_layer_mssts_f_byte_size = transfer_struct_c%ts_ls_layer_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_layer_mssts_f_shapes(1) = size(fortran_type_f%ts_ls_layer_mssts, 1)
  transfer_struct_c%ts_ls_layer_mssts_f_shapes(2) = size(fortran_type_f%ts_ls_layer_mssts, 2)
  transfer_struct_c%ts_ls_layer_mssts_f_shapes(3) = size(fortran_type_f%ts_ls_layer_mssts, 3)
  
  
  fortran_type_f%ts_ls_surf_mssts = 0_fpk
  ts_ls_surf_mssts_ptr => fortran_type_f%ts_ls_surf_mssts
  transfer_struct_c%ts_ls_surf_mssts = c_loc(ts_ls_surf_mssts_ptr(&
    lbound(fortran_type_f%ts_ls_surf_mssts,1),&
    lbound(fortran_type_f%ts_ls_surf_mssts,2)))
  inquire(iolength=transfer_struct_c%ts_ls_surf_mssts_f_byte_size) fortran_type_f%ts_ls_surf_mssts(&
    lbound(fortran_type_f%ts_ls_surf_mssts,1),&
    lbound(fortran_type_f%ts_ls_surf_mssts,2))
#ifdef ifort
  transfer_struct_c%ts_ls_surf_mssts_f_byte_size = transfer_struct_c%ts_ls_surf_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_surf_mssts_f_shapes(1) = size(fortran_type_f%ts_ls_surf_mssts, 1)
  transfer_struct_c%ts_ls_surf_mssts_f_shapes(2) = size(fortran_type_f%ts_ls_surf_mssts, 2)
  
  
end subroutine lidort_linsurf_c_init_only

subroutine lidort_linsurf_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linsurf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsurf_c_destroy

subroutine lidort_linsurf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linsurf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsurf), pointer :: fortran_type_f_from
  type(lidort_linsurf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_surfacewf = fortran_type_f_from%ts_surfacewf
  fortran_type_f_to%ts_meani_diffuse_surfwf = fortran_type_f_from%ts_meani_diffuse_surfwf
  fortran_type_f_to%ts_flux_diffuse_surfwf = fortran_type_f_from%ts_flux_diffuse_surfwf
  fortran_type_f_to%ts_sbbwfs_jacobians = fortran_type_f_from%ts_sbbwfs_jacobians
  fortran_type_f_to%ts_sbbwfs_fluxes = fortran_type_f_from%ts_sbbwfs_fluxes
  fortran_type_f_to%ts_ls_layer_mssts = fortran_type_f_from%ts_ls_layer_mssts
  fortran_type_f_to%ts_ls_surf_mssts = fortran_type_f_from%ts_ls_surf_mssts
  

end subroutine lidort_linsurf_c_copy

! Links to type: "lidort_linoutputs" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linoutputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linoutputs

  type(lidort_linoutputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linoutputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linoutputs_c_alloc_init

! Links to type: "lidort_linoutputs" from module: "lidort_lin_outputs_def_m" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linoutputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m
  use lidort_pars_m

  type(lidort_linoutputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  type(lidort_linatmos), pointer :: atmos_ptr
  type(lidort_linsurf), pointer :: surf_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  atmos_ptr => fortran_type_f%atmos
  transfer_struct_c%atmos = c_loc(atmos_ptr)
  inquire(iolength=transfer_struct_c%atmos_f_byte_size) fortran_type_f%atmos
#ifdef ifort
  transfer_struct_c%atmos_f_byte_size = transfer_struct_c%atmos_f_byte_size * 4
#endif
  
  
  surf_ptr => fortran_type_f%surf
  transfer_struct_c%surf = c_loc(surf_ptr)
  inquire(iolength=transfer_struct_c%surf_f_byte_size) fortran_type_f%surf
#ifdef ifort
  transfer_struct_c%surf_f_byte_size = transfer_struct_c%surf_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_linoutputs_c_init_only

subroutine lidort_linoutputs_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linoutputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linoutputs_c_destroy

subroutine lidort_linoutputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_outputs_def_m, only : lidort_linoutputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linoutputs), pointer :: fortran_type_f_from
  type(lidort_linoutputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%atmos = fortran_type_f_from%atmos
  fortran_type_f_to%surf = fortran_type_f_from%surf
  

end subroutine lidort_linoutputs_c_copy

! Links to type: "lidort_linsup_brdf" from module: "lidort_lin_sup_brdf_def_m" in file: "lidort_lin_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_brdf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_brdf_def_m, only : lidort_linsup_brdf

  type(lidort_linsup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_brdf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_brdf_c_alloc_init

! Links to type: "lidort_linsup_brdf" from module: "lidort_lin_sup_brdf_def_m" in file: "lidort_lin_sup_brdf_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_brdf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_brdf_def_m
  use lidort_pars_m

  type(lidort_linsup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_brdf_f_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_user_brdf_f_ptr
  real(c_double), dimension(:,:), pointer :: ts_ls_emissivity_ptr
  real(c_double), dimension(:,:), pointer :: ts_ls_user_emissivity_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_ls_exactdb_brdfunc = 0_fpk
  ts_ls_exactdb_brdfunc_ptr => fortran_type_f%ts_ls_exactdb_brdfunc
  transfer_struct_c%ts_ls_exactdb_brdfunc = c_loc(ts_ls_exactdb_brdfunc_ptr(&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,1),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,2),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,3),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,4)))
  inquire(iolength=transfer_struct_c%ts_ls_exactdb_brdfunc_f_byte_size) fortran_type_f%ts_ls_exactdb_brdfunc(&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,1),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,2),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,3),&
    lbound(fortran_type_f%ts_ls_exactdb_brdfunc,4))
#ifdef ifort
  transfer_struct_c%ts_ls_exactdb_brdfunc_f_byte_size = transfer_struct_c%ts_ls_exactdb_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_exactdb_brdfunc_f_shapes(1) = size(fortran_type_f%ts_ls_exactdb_brdfunc, 1)
  transfer_struct_c%ts_ls_exactdb_brdfunc_f_shapes(2) = size(fortran_type_f%ts_ls_exactdb_brdfunc, 2)
  transfer_struct_c%ts_ls_exactdb_brdfunc_f_shapes(3) = size(fortran_type_f%ts_ls_exactdb_brdfunc, 3)
  transfer_struct_c%ts_ls_exactdb_brdfunc_f_shapes(4) = size(fortran_type_f%ts_ls_exactdb_brdfunc, 4)
  
  
  fortran_type_f%ts_ls_brdf_f_0 = 0_fpk
  ts_ls_brdf_f_0_ptr => fortran_type_f%ts_ls_brdf_f_0
  transfer_struct_c%ts_ls_brdf_f_0 = c_loc(ts_ls_brdf_f_0_ptr(&
    lbound(fortran_type_f%ts_ls_brdf_f_0,1),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,2),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,3),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,4)))
  inquire(iolength=transfer_struct_c%ts_ls_brdf_f_0_f_byte_size) fortran_type_f%ts_ls_brdf_f_0(&
    lbound(fortran_type_f%ts_ls_brdf_f_0,1),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,2),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,3),&
    lbound(fortran_type_f%ts_ls_brdf_f_0,4))
#ifdef ifort
  transfer_struct_c%ts_ls_brdf_f_0_f_byte_size = transfer_struct_c%ts_ls_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_brdf_f_0_f_shapes(1) = size(fortran_type_f%ts_ls_brdf_f_0, 1)
  transfer_struct_c%ts_ls_brdf_f_0_f_shapes(2) = size(fortran_type_f%ts_ls_brdf_f_0, 2)
  transfer_struct_c%ts_ls_brdf_f_0_f_shapes(3) = size(fortran_type_f%ts_ls_brdf_f_0, 3)
  transfer_struct_c%ts_ls_brdf_f_0_f_shapes(4) = size(fortran_type_f%ts_ls_brdf_f_0, 4)
  
  
  fortran_type_f%ts_ls_brdf_f = 0_fpk
  ts_ls_brdf_f_ptr => fortran_type_f%ts_ls_brdf_f
  transfer_struct_c%ts_ls_brdf_f = c_loc(ts_ls_brdf_f_ptr(&
    lbound(fortran_type_f%ts_ls_brdf_f,1),&
    lbound(fortran_type_f%ts_ls_brdf_f,2),&
    lbound(fortran_type_f%ts_ls_brdf_f,3),&
    lbound(fortran_type_f%ts_ls_brdf_f,4)))
  inquire(iolength=transfer_struct_c%ts_ls_brdf_f_f_byte_size) fortran_type_f%ts_ls_brdf_f(&
    lbound(fortran_type_f%ts_ls_brdf_f,1),&
    lbound(fortran_type_f%ts_ls_brdf_f,2),&
    lbound(fortran_type_f%ts_ls_brdf_f,3),&
    lbound(fortran_type_f%ts_ls_brdf_f,4))
#ifdef ifort
  transfer_struct_c%ts_ls_brdf_f_f_byte_size = transfer_struct_c%ts_ls_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_brdf_f_f_shapes(1) = size(fortran_type_f%ts_ls_brdf_f, 1)
  transfer_struct_c%ts_ls_brdf_f_f_shapes(2) = size(fortran_type_f%ts_ls_brdf_f, 2)
  transfer_struct_c%ts_ls_brdf_f_f_shapes(3) = size(fortran_type_f%ts_ls_brdf_f, 3)
  transfer_struct_c%ts_ls_brdf_f_f_shapes(4) = size(fortran_type_f%ts_ls_brdf_f, 4)
  
  
  fortran_type_f%ts_ls_user_brdf_f_0 = 0_fpk
  ts_ls_user_brdf_f_0_ptr => fortran_type_f%ts_ls_user_brdf_f_0
  transfer_struct_c%ts_ls_user_brdf_f_0 = c_loc(ts_ls_user_brdf_f_0_ptr(&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,1),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,2),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,3),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,4)))
  inquire(iolength=transfer_struct_c%ts_ls_user_brdf_f_0_f_byte_size) fortran_type_f%ts_ls_user_brdf_f_0(&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,1),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,2),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,3),&
    lbound(fortran_type_f%ts_ls_user_brdf_f_0,4))
#ifdef ifort
  transfer_struct_c%ts_ls_user_brdf_f_0_f_byte_size = transfer_struct_c%ts_ls_user_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_user_brdf_f_0_f_shapes(1) = size(fortran_type_f%ts_ls_user_brdf_f_0, 1)
  transfer_struct_c%ts_ls_user_brdf_f_0_f_shapes(2) = size(fortran_type_f%ts_ls_user_brdf_f_0, 2)
  transfer_struct_c%ts_ls_user_brdf_f_0_f_shapes(3) = size(fortran_type_f%ts_ls_user_brdf_f_0, 3)
  transfer_struct_c%ts_ls_user_brdf_f_0_f_shapes(4) = size(fortran_type_f%ts_ls_user_brdf_f_0, 4)
  
  
  fortran_type_f%ts_ls_user_brdf_f = 0_fpk
  ts_ls_user_brdf_f_ptr => fortran_type_f%ts_ls_user_brdf_f
  transfer_struct_c%ts_ls_user_brdf_f = c_loc(ts_ls_user_brdf_f_ptr(&
    lbound(fortran_type_f%ts_ls_user_brdf_f,1),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,2),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,3),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,4)))
  inquire(iolength=transfer_struct_c%ts_ls_user_brdf_f_f_byte_size) fortran_type_f%ts_ls_user_brdf_f(&
    lbound(fortran_type_f%ts_ls_user_brdf_f,1),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,2),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,3),&
    lbound(fortran_type_f%ts_ls_user_brdf_f,4))
#ifdef ifort
  transfer_struct_c%ts_ls_user_brdf_f_f_byte_size = transfer_struct_c%ts_ls_user_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_user_brdf_f_f_shapes(1) = size(fortran_type_f%ts_ls_user_brdf_f, 1)
  transfer_struct_c%ts_ls_user_brdf_f_f_shapes(2) = size(fortran_type_f%ts_ls_user_brdf_f, 2)
  transfer_struct_c%ts_ls_user_brdf_f_f_shapes(3) = size(fortran_type_f%ts_ls_user_brdf_f, 3)
  transfer_struct_c%ts_ls_user_brdf_f_f_shapes(4) = size(fortran_type_f%ts_ls_user_brdf_f, 4)
  
  
  fortran_type_f%ts_ls_emissivity = 0_fpk
  ts_ls_emissivity_ptr => fortran_type_f%ts_ls_emissivity
  transfer_struct_c%ts_ls_emissivity = c_loc(ts_ls_emissivity_ptr(&
    lbound(fortran_type_f%ts_ls_emissivity,1),&
    lbound(fortran_type_f%ts_ls_emissivity,2)))
  inquire(iolength=transfer_struct_c%ts_ls_emissivity_f_byte_size) fortran_type_f%ts_ls_emissivity(&
    lbound(fortran_type_f%ts_ls_emissivity,1),&
    lbound(fortran_type_f%ts_ls_emissivity,2))
#ifdef ifort
  transfer_struct_c%ts_ls_emissivity_f_byte_size = transfer_struct_c%ts_ls_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_emissivity_f_shapes(1) = size(fortran_type_f%ts_ls_emissivity, 1)
  transfer_struct_c%ts_ls_emissivity_f_shapes(2) = size(fortran_type_f%ts_ls_emissivity, 2)
  
  
  fortran_type_f%ts_ls_user_emissivity = 0_fpk
  ts_ls_user_emissivity_ptr => fortran_type_f%ts_ls_user_emissivity
  transfer_struct_c%ts_ls_user_emissivity = c_loc(ts_ls_user_emissivity_ptr(&
    lbound(fortran_type_f%ts_ls_user_emissivity,1),&
    lbound(fortran_type_f%ts_ls_user_emissivity,2)))
  inquire(iolength=transfer_struct_c%ts_ls_user_emissivity_f_byte_size) fortran_type_f%ts_ls_user_emissivity(&
    lbound(fortran_type_f%ts_ls_user_emissivity,1),&
    lbound(fortran_type_f%ts_ls_user_emissivity,2))
#ifdef ifort
  transfer_struct_c%ts_ls_user_emissivity_f_byte_size = transfer_struct_c%ts_ls_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_user_emissivity_f_shapes(1) = size(fortran_type_f%ts_ls_user_emissivity, 1)
  transfer_struct_c%ts_ls_user_emissivity_f_shapes(2) = size(fortran_type_f%ts_ls_user_emissivity, 2)
  
  
end subroutine lidort_linsup_brdf_c_init_only

subroutine lidort_linsup_brdf_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_brdf_def_m, only : lidort_linsup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_brdf_c_destroy

subroutine lidort_linsup_brdf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_brdf_def_m, only : lidort_linsup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_brdf), pointer :: fortran_type_f_from
  type(lidort_linsup_brdf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_ls_exactdb_brdfunc = fortran_type_f_from%ts_ls_exactdb_brdfunc
  fortran_type_f_to%ts_ls_brdf_f_0 = fortran_type_f_from%ts_ls_brdf_f_0
  fortran_type_f_to%ts_ls_brdf_f = fortran_type_f_from%ts_ls_brdf_f
  fortran_type_f_to%ts_ls_user_brdf_f_0 = fortran_type_f_from%ts_ls_user_brdf_f_0
  fortran_type_f_to%ts_ls_user_brdf_f = fortran_type_f_from%ts_ls_user_brdf_f
  fortran_type_f_to%ts_ls_emissivity = fortran_type_f_from%ts_ls_emissivity
  fortran_type_f_to%ts_ls_user_emissivity = fortran_type_f_from%ts_ls_user_emissivity
  

end subroutine lidort_linsup_brdf_c_copy

! Links to type: "lidort_linsup_sleave" from module: "lidort_lin_sup_sleave_def_m" in file: "lidort_lin_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_sleave_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_sleave_def_m, only : lidort_linsup_sleave

  type(lidort_linsup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_sleave_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_sleave_c_alloc_init

! Links to type: "lidort_linsup_sleave" from module: "lidort_lin_sup_sleave_def_m" in file: "lidort_lin_sup_sleave_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_sleave_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_sleave_def_m
  use lidort_pars_m

  type(lidort_linsup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  real(c_double), dimension(:,:), pointer :: ts_lssl_slterm_isotropic_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_lssl_slterm_userangles_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_lssl_slterm_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_lssl_user_slterm_f_0_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_lssl_slterm_isotropic = 0_fpk
  ts_lssl_slterm_isotropic_ptr => fortran_type_f%ts_lssl_slterm_isotropic
  transfer_struct_c%ts_lssl_slterm_isotropic = c_loc(ts_lssl_slterm_isotropic_ptr(&
    lbound(fortran_type_f%ts_lssl_slterm_isotropic,1),&
    lbound(fortran_type_f%ts_lssl_slterm_isotropic,2)))
  inquire(iolength=transfer_struct_c%ts_lssl_slterm_isotropic_f_byte_size) fortran_type_f%ts_lssl_slterm_isotropic(&
    lbound(fortran_type_f%ts_lssl_slterm_isotropic,1),&
    lbound(fortran_type_f%ts_lssl_slterm_isotropic,2))
#ifdef ifort
  transfer_struct_c%ts_lssl_slterm_isotropic_f_byte_size = transfer_struct_c%ts_lssl_slterm_isotropic_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lssl_slterm_isotropic_f_shapes(1) = size(fortran_type_f%ts_lssl_slterm_isotropic, 1)
  transfer_struct_c%ts_lssl_slterm_isotropic_f_shapes(2) = size(fortran_type_f%ts_lssl_slterm_isotropic, 2)
  
  
  fortran_type_f%ts_lssl_slterm_userangles = 0_fpk
  ts_lssl_slterm_userangles_ptr => fortran_type_f%ts_lssl_slterm_userangles
  transfer_struct_c%ts_lssl_slterm_userangles = c_loc(ts_lssl_slterm_userangles_ptr(&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,1),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,2),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,3),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,4)))
  inquire(iolength=transfer_struct_c%ts_lssl_slterm_userangles_f_byte_size) fortran_type_f%ts_lssl_slterm_userangles(&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,1),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,2),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,3),&
    lbound(fortran_type_f%ts_lssl_slterm_userangles,4))
#ifdef ifort
  transfer_struct_c%ts_lssl_slterm_userangles_f_byte_size = transfer_struct_c%ts_lssl_slterm_userangles_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lssl_slterm_userangles_f_shapes(1) = size(fortran_type_f%ts_lssl_slterm_userangles, 1)
  transfer_struct_c%ts_lssl_slterm_userangles_f_shapes(2) = size(fortran_type_f%ts_lssl_slterm_userangles, 2)
  transfer_struct_c%ts_lssl_slterm_userangles_f_shapes(3) = size(fortran_type_f%ts_lssl_slterm_userangles, 3)
  transfer_struct_c%ts_lssl_slterm_userangles_f_shapes(4) = size(fortran_type_f%ts_lssl_slterm_userangles, 4)
  
  
  fortran_type_f%ts_lssl_slterm_f_0 = 0_fpk
  ts_lssl_slterm_f_0_ptr => fortran_type_f%ts_lssl_slterm_f_0
  transfer_struct_c%ts_lssl_slterm_f_0 = c_loc(ts_lssl_slterm_f_0_ptr(&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,1),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,2),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,3),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,4)))
  inquire(iolength=transfer_struct_c%ts_lssl_slterm_f_0_f_byte_size) fortran_type_f%ts_lssl_slterm_f_0(&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,1),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,2),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,3),&
    lbound(fortran_type_f%ts_lssl_slterm_f_0,4))
#ifdef ifort
  transfer_struct_c%ts_lssl_slterm_f_0_f_byte_size = transfer_struct_c%ts_lssl_slterm_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lssl_slterm_f_0_f_shapes(1) = size(fortran_type_f%ts_lssl_slterm_f_0, 1)
  transfer_struct_c%ts_lssl_slterm_f_0_f_shapes(2) = size(fortran_type_f%ts_lssl_slterm_f_0, 2)
  transfer_struct_c%ts_lssl_slterm_f_0_f_shapes(3) = size(fortran_type_f%ts_lssl_slterm_f_0, 3)
  transfer_struct_c%ts_lssl_slterm_f_0_f_shapes(4) = size(fortran_type_f%ts_lssl_slterm_f_0, 4)
  
  
  fortran_type_f%ts_lssl_user_slterm_f_0 = 0_fpk
  ts_lssl_user_slterm_f_0_ptr => fortran_type_f%ts_lssl_user_slterm_f_0
  transfer_struct_c%ts_lssl_user_slterm_f_0 = c_loc(ts_lssl_user_slterm_f_0_ptr(&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,1),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,2),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,3),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,4)))
  inquire(iolength=transfer_struct_c%ts_lssl_user_slterm_f_0_f_byte_size) fortran_type_f%ts_lssl_user_slterm_f_0(&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,1),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,2),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,3),&
    lbound(fortran_type_f%ts_lssl_user_slterm_f_0,4))
#ifdef ifort
  transfer_struct_c%ts_lssl_user_slterm_f_0_f_byte_size = transfer_struct_c%ts_lssl_user_slterm_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lssl_user_slterm_f_0_f_shapes(1) = size(fortran_type_f%ts_lssl_user_slterm_f_0, 1)
  transfer_struct_c%ts_lssl_user_slterm_f_0_f_shapes(2) = size(fortran_type_f%ts_lssl_user_slterm_f_0, 2)
  transfer_struct_c%ts_lssl_user_slterm_f_0_f_shapes(3) = size(fortran_type_f%ts_lssl_user_slterm_f_0, 3)
  transfer_struct_c%ts_lssl_user_slterm_f_0_f_shapes(4) = size(fortran_type_f%ts_lssl_user_slterm_f_0, 4)
  
  
end subroutine lidort_linsup_sleave_c_init_only

subroutine lidort_linsup_sleave_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_sleave_def_m, only : lidort_linsup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_sleave_c_destroy

subroutine lidort_linsup_sleave_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_sleave_def_m, only : lidort_linsup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_sleave), pointer :: fortran_type_f_from
  type(lidort_linsup_sleave), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_lssl_slterm_isotropic = fortran_type_f_from%ts_lssl_slterm_isotropic
  fortran_type_f_to%ts_lssl_slterm_userangles = fortran_type_f_from%ts_lssl_slterm_userangles
  fortran_type_f_to%ts_lssl_slterm_f_0 = fortran_type_f_from%ts_lssl_slterm_f_0
  fortran_type_f_to%ts_lssl_user_slterm_f_0 = fortran_type_f_from%ts_lssl_user_slterm_f_0
  

end subroutine lidort_linsup_sleave_c_copy

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_atmos_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_atmos

  type(lidort_linsup_ss_atmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_atmos_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_atmos_c_alloc_init

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_atmos_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m
  use lidort_pars_m

  type(lidort_linsup_ss_atmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_columnwf_ss_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_columnwf_db_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_profilewf_ss_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_profilewf_db_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_columnwf_ss = 0_fpk
  ts_columnwf_ss_ptr => fortran_type_f%ts_columnwf_ss
  transfer_struct_c%ts_columnwf_ss = c_loc(ts_columnwf_ss_ptr(&
    lbound(fortran_type_f%ts_columnwf_ss,1),&
    lbound(fortran_type_f%ts_columnwf_ss,2),&
    lbound(fortran_type_f%ts_columnwf_ss,3),&
    lbound(fortran_type_f%ts_columnwf_ss,4)))
  inquire(iolength=transfer_struct_c%ts_columnwf_ss_f_byte_size) fortran_type_f%ts_columnwf_ss(&
    lbound(fortran_type_f%ts_columnwf_ss,1),&
    lbound(fortran_type_f%ts_columnwf_ss,2),&
    lbound(fortran_type_f%ts_columnwf_ss,3),&
    lbound(fortran_type_f%ts_columnwf_ss,4))
#ifdef ifort
  transfer_struct_c%ts_columnwf_ss_f_byte_size = transfer_struct_c%ts_columnwf_ss_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_columnwf_ss_f_shapes(1) = size(fortran_type_f%ts_columnwf_ss, 1)
  transfer_struct_c%ts_columnwf_ss_f_shapes(2) = size(fortran_type_f%ts_columnwf_ss, 2)
  transfer_struct_c%ts_columnwf_ss_f_shapes(3) = size(fortran_type_f%ts_columnwf_ss, 3)
  transfer_struct_c%ts_columnwf_ss_f_shapes(4) = size(fortran_type_f%ts_columnwf_ss, 4)
  
  
  fortran_type_f%ts_columnwf_db = 0_fpk
  ts_columnwf_db_ptr => fortran_type_f%ts_columnwf_db
  transfer_struct_c%ts_columnwf_db = c_loc(ts_columnwf_db_ptr(&
    lbound(fortran_type_f%ts_columnwf_db,1),&
    lbound(fortran_type_f%ts_columnwf_db,2),&
    lbound(fortran_type_f%ts_columnwf_db,3)))
  inquire(iolength=transfer_struct_c%ts_columnwf_db_f_byte_size) fortran_type_f%ts_columnwf_db(&
    lbound(fortran_type_f%ts_columnwf_db,1),&
    lbound(fortran_type_f%ts_columnwf_db,2),&
    lbound(fortran_type_f%ts_columnwf_db,3))
#ifdef ifort
  transfer_struct_c%ts_columnwf_db_f_byte_size = transfer_struct_c%ts_columnwf_db_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_columnwf_db_f_shapes(1) = size(fortran_type_f%ts_columnwf_db, 1)
  transfer_struct_c%ts_columnwf_db_f_shapes(2) = size(fortran_type_f%ts_columnwf_db, 2)
  transfer_struct_c%ts_columnwf_db_f_shapes(3) = size(fortran_type_f%ts_columnwf_db, 3)
  
  
  fortran_type_f%ts_profilewf_ss = 0_fpk
  ts_profilewf_ss_ptr => fortran_type_f%ts_profilewf_ss
  transfer_struct_c%ts_profilewf_ss = c_loc(ts_profilewf_ss_ptr(&
    lbound(fortran_type_f%ts_profilewf_ss,1),&
    lbound(fortran_type_f%ts_profilewf_ss,2),&
    lbound(fortran_type_f%ts_profilewf_ss,3),&
    lbound(fortran_type_f%ts_profilewf_ss,4),&
    lbound(fortran_type_f%ts_profilewf_ss,5)))
  inquire(iolength=transfer_struct_c%ts_profilewf_ss_f_byte_size) fortran_type_f%ts_profilewf_ss(&
    lbound(fortran_type_f%ts_profilewf_ss,1),&
    lbound(fortran_type_f%ts_profilewf_ss,2),&
    lbound(fortran_type_f%ts_profilewf_ss,3),&
    lbound(fortran_type_f%ts_profilewf_ss,4),&
    lbound(fortran_type_f%ts_profilewf_ss,5))
#ifdef ifort
  transfer_struct_c%ts_profilewf_ss_f_byte_size = transfer_struct_c%ts_profilewf_ss_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_profilewf_ss_f_shapes(1) = size(fortran_type_f%ts_profilewf_ss, 1)
  transfer_struct_c%ts_profilewf_ss_f_shapes(2) = size(fortran_type_f%ts_profilewf_ss, 2)
  transfer_struct_c%ts_profilewf_ss_f_shapes(3) = size(fortran_type_f%ts_profilewf_ss, 3)
  transfer_struct_c%ts_profilewf_ss_f_shapes(4) = size(fortran_type_f%ts_profilewf_ss, 4)
  transfer_struct_c%ts_profilewf_ss_f_shapes(5) = size(fortran_type_f%ts_profilewf_ss, 5)
  
  
  fortran_type_f%ts_profilewf_db = 0_fpk
  ts_profilewf_db_ptr => fortran_type_f%ts_profilewf_db
  transfer_struct_c%ts_profilewf_db = c_loc(ts_profilewf_db_ptr(&
    lbound(fortran_type_f%ts_profilewf_db,1),&
    lbound(fortran_type_f%ts_profilewf_db,2),&
    lbound(fortran_type_f%ts_profilewf_db,3),&
    lbound(fortran_type_f%ts_profilewf_db,4)))
  inquire(iolength=transfer_struct_c%ts_profilewf_db_f_byte_size) fortran_type_f%ts_profilewf_db(&
    lbound(fortran_type_f%ts_profilewf_db,1),&
    lbound(fortran_type_f%ts_profilewf_db,2),&
    lbound(fortran_type_f%ts_profilewf_db,3),&
    lbound(fortran_type_f%ts_profilewf_db,4))
#ifdef ifort
  transfer_struct_c%ts_profilewf_db_f_byte_size = transfer_struct_c%ts_profilewf_db_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_profilewf_db_f_shapes(1) = size(fortran_type_f%ts_profilewf_db, 1)
  transfer_struct_c%ts_profilewf_db_f_shapes(2) = size(fortran_type_f%ts_profilewf_db, 2)
  transfer_struct_c%ts_profilewf_db_f_shapes(3) = size(fortran_type_f%ts_profilewf_db, 3)
  transfer_struct_c%ts_profilewf_db_f_shapes(4) = size(fortran_type_f%ts_profilewf_db, 4)
  
  
end subroutine lidort_linsup_ss_atmos_c_init_only

subroutine lidort_linsup_ss_atmos_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_atmos

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_atmos_c_destroy

subroutine lidort_linsup_ss_atmos_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_atmos

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f_from
  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_columnwf_ss = fortran_type_f_from%ts_columnwf_ss
  fortran_type_f_to%ts_columnwf_db = fortran_type_f_from%ts_columnwf_db
  fortran_type_f_to%ts_profilewf_ss = fortran_type_f_from%ts_profilewf_ss
  fortran_type_f_to%ts_profilewf_db = fortran_type_f_from%ts_profilewf_db
  

end subroutine lidort_linsup_ss_atmos_c_copy

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_surf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_surf

  type(lidort_linsup_ss_surf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_surf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_surf_c_alloc_init

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_surf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m
  use lidort_pars_m

  type(lidort_linsup_ss_surf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_surfacewf_db_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_surfacewf_db = 0_fpk
  ts_surfacewf_db_ptr => fortran_type_f%ts_surfacewf_db
  transfer_struct_c%ts_surfacewf_db = c_loc(ts_surfacewf_db_ptr(&
    lbound(fortran_type_f%ts_surfacewf_db,1),&
    lbound(fortran_type_f%ts_surfacewf_db,2),&
    lbound(fortran_type_f%ts_surfacewf_db,3)))
  inquire(iolength=transfer_struct_c%ts_surfacewf_db_f_byte_size) fortran_type_f%ts_surfacewf_db(&
    lbound(fortran_type_f%ts_surfacewf_db,1),&
    lbound(fortran_type_f%ts_surfacewf_db,2),&
    lbound(fortran_type_f%ts_surfacewf_db,3))
#ifdef ifort
  transfer_struct_c%ts_surfacewf_db_f_byte_size = transfer_struct_c%ts_surfacewf_db_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_surfacewf_db_f_shapes(1) = size(fortran_type_f%ts_surfacewf_db, 1)
  transfer_struct_c%ts_surfacewf_db_f_shapes(2) = size(fortran_type_f%ts_surfacewf_db, 2)
  transfer_struct_c%ts_surfacewf_db_f_shapes(3) = size(fortran_type_f%ts_surfacewf_db, 3)
  
  
end subroutine lidort_linsup_ss_surf_c_init_only

subroutine lidort_linsup_ss_surf_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_surf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_surf_c_destroy

subroutine lidort_linsup_ss_surf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss_surf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f_from
  type(lidort_linsup_ss_surf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_surfacewf_db = fortran_type_f_from%ts_surfacewf_db
  

end subroutine lidort_linsup_ss_surf_c_copy

! Links to type: "lidort_linsup_ss" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss

  type(lidort_linsup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_c_alloc_init

! Links to type: "lidort_linsup_ss" from module: "lidort_lin_sup_ss_def_m" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m
  use lidort_pars_m

  type(lidort_linsup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  type(lidort_linsup_ss_atmos), pointer :: atmos_ptr
  type(lidort_linsup_ss_surf), pointer :: surf_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  atmos_ptr => fortran_type_f%atmos
  transfer_struct_c%atmos = c_loc(atmos_ptr)
  inquire(iolength=transfer_struct_c%atmos_f_byte_size) fortran_type_f%atmos
#ifdef ifort
  transfer_struct_c%atmos_f_byte_size = transfer_struct_c%atmos_f_byte_size * 4
#endif
  
  
  surf_ptr => fortran_type_f%surf
  transfer_struct_c%surf = c_loc(surf_ptr)
  inquire(iolength=transfer_struct_c%surf_f_byte_size) fortran_type_f%surf
#ifdef ifort
  transfer_struct_c%surf_f_byte_size = transfer_struct_c%surf_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_linsup_ss_c_init_only

subroutine lidort_linsup_ss_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_c_destroy

subroutine lidort_linsup_ss_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_ss_def_m, only : lidort_linsup_ss

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_ss), pointer :: fortran_type_f_from
  type(lidort_linsup_ss), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%atmos = fortran_type_f_from%atmos
  fortran_type_f_to%surf = fortran_type_f_from%surf
  

end subroutine lidort_linsup_ss_c_copy

! Links to type: "lidort_linsup_inout" from module: "lidort_lin_sup_inout_def_m" in file: "lidort_lin_sup_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_inout_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_inout_def_m, only : lidort_linsup_inout

  type(lidort_linsup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_inout_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_inout_c_alloc_init

! Links to type: "lidort_linsup_inout" from module: "lidort_lin_sup_inout_def_m" in file: "lidort_lin_sup_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_inout_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lin_sup_inout_def_m
  use lidort_lin_sup_brdf_def_m
  use lidort_lin_sup_ss_def_m
  use lidort_lin_sup_sleave_def_m

  type(lidort_linsup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  type(lidort_linsup_brdf), pointer :: brdf_ptr
  type(lidort_linsup_ss), pointer :: ss_ptr
  type(lidort_linsup_sleave), pointer :: sleave_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  brdf_ptr => fortran_type_f%brdf
  transfer_struct_c%brdf = c_loc(brdf_ptr)
  inquire(iolength=transfer_struct_c%brdf_f_byte_size) fortran_type_f%brdf
#ifdef ifort
  transfer_struct_c%brdf_f_byte_size = transfer_struct_c%brdf_f_byte_size * 4
#endif
  
  
  ss_ptr => fortran_type_f%ss
  transfer_struct_c%ss = c_loc(ss_ptr)
  inquire(iolength=transfer_struct_c%ss_f_byte_size) fortran_type_f%ss
#ifdef ifort
  transfer_struct_c%ss_f_byte_size = transfer_struct_c%ss_f_byte_size * 4
#endif
  
  
  sleave_ptr => fortran_type_f%sleave
  transfer_struct_c%sleave = c_loc(sleave_ptr)
  inquire(iolength=transfer_struct_c%sleave_f_byte_size) fortran_type_f%sleave
#ifdef ifort
  transfer_struct_c%sleave_f_byte_size = transfer_struct_c%sleave_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_linsup_inout_c_init_only

subroutine lidort_linsup_inout_c_destroy(fortran_type_c) bind(C)
  use lidort_lin_sup_inout_def_m, only : lidort_linsup_inout

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_inout_c_destroy

subroutine lidort_linsup_inout_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lin_sup_inout_def_m, only : lidort_linsup_inout

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_inout), pointer :: fortran_type_f_from
  type(lidort_linsup_inout), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%brdf = fortran_type_f_from%brdf
  fortran_type_f_to%ss = fortran_type_f_from%ss
  fortran_type_f_to%sleave = fortran_type_f_from%sleave
  

end subroutine lidort_linsup_inout_c_copy

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_main_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_main_outputs

  type(lidort_main_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_main_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_main_outputs_c_alloc_init

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_main_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m
  use lidort_pars_m

  type(lidort_main_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_intensity_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_meani_diffuse_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_flux_diffuse_ptr
  real(c_double), dimension(:,:), pointer :: ts_dnmeani_direct_ptr
  real(c_double), dimension(:,:), pointer :: ts_dnflux_direct_ptr
  real(c_double), dimension(:), pointer :: ts_albmed_user_ptr
  real(c_double), dimension(:), pointer :: ts_trnmed_user_ptr
  real(c_double), dimension(:), pointer :: ts_albmed_fluxes_ptr
  real(c_double), dimension(:), pointer :: ts_trnmed_fluxes_ptr
  real(c_double), dimension(:), pointer :: ts_planetary_transterm_ptr
  real(c_double), pointer :: ts_planetary_sbterm_ptr
  real(c_double), dimension(:,:), pointer :: ts_pathgeoms_ptr
  real(c_double), dimension(:,:), pointer :: ts_lostrans_ptr
  real(c_double), dimension(:,:), pointer :: ts_layer_mssts_ptr
  real(c_double), dimension(:), pointer :: ts_surf_mssts_ptr
  real(c_double), dimension(:,:), pointer :: ts_contribs_ptr
  integer(c_int), dimension(:), pointer :: ts_fourier_saved_ptr
  integer(c_int), pointer :: ts_n_geometries_ptr
  real(c_double), dimension(:), pointer :: ts_solarbeam_boatrans_ptr
  real(c_double), pointer :: ts_spheralb_ptr
  real(c_double), dimension(:), pointer :: ts_trans1_user_ptr
  real(c_double), dimension(:), pointer :: ts_trans1_beam_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_intensity = 0_fpk
  ts_intensity_ptr => fortran_type_f%ts_intensity
  transfer_struct_c%ts_intensity = c_loc(ts_intensity_ptr(&
    lbound(fortran_type_f%ts_intensity,1),&
    lbound(fortran_type_f%ts_intensity,2),&
    lbound(fortran_type_f%ts_intensity,3)))
  inquire(iolength=transfer_struct_c%ts_intensity_f_byte_size) fortran_type_f%ts_intensity(&
    lbound(fortran_type_f%ts_intensity,1),&
    lbound(fortran_type_f%ts_intensity,2),&
    lbound(fortran_type_f%ts_intensity,3))
#ifdef ifort
  transfer_struct_c%ts_intensity_f_byte_size = transfer_struct_c%ts_intensity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_intensity_f_shapes(1) = size(fortran_type_f%ts_intensity, 1)
  transfer_struct_c%ts_intensity_f_shapes(2) = size(fortran_type_f%ts_intensity, 2)
  transfer_struct_c%ts_intensity_f_shapes(3) = size(fortran_type_f%ts_intensity, 3)
  
  
  fortran_type_f%ts_meani_diffuse = 0_fpk
  ts_meani_diffuse_ptr => fortran_type_f%ts_meani_diffuse
  transfer_struct_c%ts_meani_diffuse = c_loc(ts_meani_diffuse_ptr(&
    lbound(fortran_type_f%ts_meani_diffuse,1),&
    lbound(fortran_type_f%ts_meani_diffuse,2),&
    lbound(fortran_type_f%ts_meani_diffuse,3)))
  inquire(iolength=transfer_struct_c%ts_meani_diffuse_f_byte_size) fortran_type_f%ts_meani_diffuse(&
    lbound(fortran_type_f%ts_meani_diffuse,1),&
    lbound(fortran_type_f%ts_meani_diffuse,2),&
    lbound(fortran_type_f%ts_meani_diffuse,3))
#ifdef ifort
  transfer_struct_c%ts_meani_diffuse_f_byte_size = transfer_struct_c%ts_meani_diffuse_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_meani_diffuse_f_shapes(1) = size(fortran_type_f%ts_meani_diffuse, 1)
  transfer_struct_c%ts_meani_diffuse_f_shapes(2) = size(fortran_type_f%ts_meani_diffuse, 2)
  transfer_struct_c%ts_meani_diffuse_f_shapes(3) = size(fortran_type_f%ts_meani_diffuse, 3)
  
  
  fortran_type_f%ts_flux_diffuse = 0_fpk
  ts_flux_diffuse_ptr => fortran_type_f%ts_flux_diffuse
  transfer_struct_c%ts_flux_diffuse = c_loc(ts_flux_diffuse_ptr(&
    lbound(fortran_type_f%ts_flux_diffuse,1),&
    lbound(fortran_type_f%ts_flux_diffuse,2),&
    lbound(fortran_type_f%ts_flux_diffuse,3)))
  inquire(iolength=transfer_struct_c%ts_flux_diffuse_f_byte_size) fortran_type_f%ts_flux_diffuse(&
    lbound(fortran_type_f%ts_flux_diffuse,1),&
    lbound(fortran_type_f%ts_flux_diffuse,2),&
    lbound(fortran_type_f%ts_flux_diffuse,3))
#ifdef ifort
  transfer_struct_c%ts_flux_diffuse_f_byte_size = transfer_struct_c%ts_flux_diffuse_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_diffuse_f_shapes(1) = size(fortran_type_f%ts_flux_diffuse, 1)
  transfer_struct_c%ts_flux_diffuse_f_shapes(2) = size(fortran_type_f%ts_flux_diffuse, 2)
  transfer_struct_c%ts_flux_diffuse_f_shapes(3) = size(fortran_type_f%ts_flux_diffuse, 3)
  
  
  fortran_type_f%ts_dnmeani_direct = 0_fpk
  ts_dnmeani_direct_ptr => fortran_type_f%ts_dnmeani_direct
  transfer_struct_c%ts_dnmeani_direct = c_loc(ts_dnmeani_direct_ptr(&
    lbound(fortran_type_f%ts_dnmeani_direct,1),&
    lbound(fortran_type_f%ts_dnmeani_direct,2)))
  inquire(iolength=transfer_struct_c%ts_dnmeani_direct_f_byte_size) fortran_type_f%ts_dnmeani_direct(&
    lbound(fortran_type_f%ts_dnmeani_direct,1),&
    lbound(fortran_type_f%ts_dnmeani_direct,2))
#ifdef ifort
  transfer_struct_c%ts_dnmeani_direct_f_byte_size = transfer_struct_c%ts_dnmeani_direct_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnmeani_direct_f_shapes(1) = size(fortran_type_f%ts_dnmeani_direct, 1)
  transfer_struct_c%ts_dnmeani_direct_f_shapes(2) = size(fortran_type_f%ts_dnmeani_direct, 2)
  
  
  fortran_type_f%ts_dnflux_direct = 0_fpk
  ts_dnflux_direct_ptr => fortran_type_f%ts_dnflux_direct
  transfer_struct_c%ts_dnflux_direct = c_loc(ts_dnflux_direct_ptr(&
    lbound(fortran_type_f%ts_dnflux_direct,1),&
    lbound(fortran_type_f%ts_dnflux_direct,2)))
  inquire(iolength=transfer_struct_c%ts_dnflux_direct_f_byte_size) fortran_type_f%ts_dnflux_direct(&
    lbound(fortran_type_f%ts_dnflux_direct,1),&
    lbound(fortran_type_f%ts_dnflux_direct,2))
#ifdef ifort
  transfer_struct_c%ts_dnflux_direct_f_byte_size = transfer_struct_c%ts_dnflux_direct_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnflux_direct_f_shapes(1) = size(fortran_type_f%ts_dnflux_direct, 1)
  transfer_struct_c%ts_dnflux_direct_f_shapes(2) = size(fortran_type_f%ts_dnflux_direct, 2)
  
  
  fortran_type_f%ts_albmed_user = 0_fpk
  ts_albmed_user_ptr => fortran_type_f%ts_albmed_user
  transfer_struct_c%ts_albmed_user = c_loc(ts_albmed_user_ptr(&
    lbound(fortran_type_f%ts_albmed_user,1)))
  inquire(iolength=transfer_struct_c%ts_albmed_user_f_byte_size) fortran_type_f%ts_albmed_user(&
    lbound(fortran_type_f%ts_albmed_user,1))
#ifdef ifort
  transfer_struct_c%ts_albmed_user_f_byte_size = transfer_struct_c%ts_albmed_user_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_user_f_shapes(1) = size(fortran_type_f%ts_albmed_user, 1)
  
  
  fortran_type_f%ts_trnmed_user = 0_fpk
  ts_trnmed_user_ptr => fortran_type_f%ts_trnmed_user
  transfer_struct_c%ts_trnmed_user = c_loc(ts_trnmed_user_ptr(&
    lbound(fortran_type_f%ts_trnmed_user,1)))
  inquire(iolength=transfer_struct_c%ts_trnmed_user_f_byte_size) fortran_type_f%ts_trnmed_user(&
    lbound(fortran_type_f%ts_trnmed_user,1))
#ifdef ifort
  transfer_struct_c%ts_trnmed_user_f_byte_size = transfer_struct_c%ts_trnmed_user_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_user_f_shapes(1) = size(fortran_type_f%ts_trnmed_user, 1)
  
  
  fortran_type_f%ts_albmed_fluxes = 0_fpk
  ts_albmed_fluxes_ptr => fortran_type_f%ts_albmed_fluxes
  transfer_struct_c%ts_albmed_fluxes = c_loc(ts_albmed_fluxes_ptr(&
    lbound(fortran_type_f%ts_albmed_fluxes,1)))
  inquire(iolength=transfer_struct_c%ts_albmed_fluxes_f_byte_size) fortran_type_f%ts_albmed_fluxes(&
    lbound(fortran_type_f%ts_albmed_fluxes,1))
#ifdef ifort
  transfer_struct_c%ts_albmed_fluxes_f_byte_size = transfer_struct_c%ts_albmed_fluxes_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_albmed_fluxes_f_shapes(1) = size(fortran_type_f%ts_albmed_fluxes, 1)
  
  
  fortran_type_f%ts_trnmed_fluxes = 0_fpk
  ts_trnmed_fluxes_ptr => fortran_type_f%ts_trnmed_fluxes
  transfer_struct_c%ts_trnmed_fluxes = c_loc(ts_trnmed_fluxes_ptr(&
    lbound(fortran_type_f%ts_trnmed_fluxes,1)))
  inquire(iolength=transfer_struct_c%ts_trnmed_fluxes_f_byte_size) fortran_type_f%ts_trnmed_fluxes(&
    lbound(fortran_type_f%ts_trnmed_fluxes,1))
#ifdef ifort
  transfer_struct_c%ts_trnmed_fluxes_f_byte_size = transfer_struct_c%ts_trnmed_fluxes_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trnmed_fluxes_f_shapes(1) = size(fortran_type_f%ts_trnmed_fluxes, 1)
  
  
  fortran_type_f%ts_planetary_transterm = 0_fpk
  ts_planetary_transterm_ptr => fortran_type_f%ts_planetary_transterm
  transfer_struct_c%ts_planetary_transterm = c_loc(ts_planetary_transterm_ptr(&
    lbound(fortran_type_f%ts_planetary_transterm,1)))
  inquire(iolength=transfer_struct_c%ts_planetary_transterm_f_byte_size) fortran_type_f%ts_planetary_transterm(&
    lbound(fortran_type_f%ts_planetary_transterm,1))
#ifdef ifort
  transfer_struct_c%ts_planetary_transterm_f_byte_size = transfer_struct_c%ts_planetary_transterm_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_planetary_transterm_f_shapes(1) = size(fortran_type_f%ts_planetary_transterm, 1)
  
  
  fortran_type_f%ts_planetary_sbterm = 0_fpk
  ts_planetary_sbterm_ptr => fortran_type_f%ts_planetary_sbterm
  transfer_struct_c%ts_planetary_sbterm = c_loc(ts_planetary_sbterm_ptr)
  inquire(iolength=transfer_struct_c%ts_planetary_sbterm_f_byte_size) fortran_type_f%ts_planetary_sbterm
#ifdef ifort
  transfer_struct_c%ts_planetary_sbterm_f_byte_size = transfer_struct_c%ts_planetary_sbterm_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_pathgeoms = 0_fpk
  ts_pathgeoms_ptr => fortran_type_f%ts_pathgeoms
  transfer_struct_c%ts_pathgeoms = c_loc(ts_pathgeoms_ptr(&
    lbound(fortran_type_f%ts_pathgeoms,1),&
    lbound(fortran_type_f%ts_pathgeoms,2)))
  inquire(iolength=transfer_struct_c%ts_pathgeoms_f_byte_size) fortran_type_f%ts_pathgeoms(&
    lbound(fortran_type_f%ts_pathgeoms,1),&
    lbound(fortran_type_f%ts_pathgeoms,2))
#ifdef ifort
  transfer_struct_c%ts_pathgeoms_f_byte_size = transfer_struct_c%ts_pathgeoms_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_pathgeoms_f_shapes(1) = size(fortran_type_f%ts_pathgeoms, 1)
  transfer_struct_c%ts_pathgeoms_f_shapes(2) = size(fortran_type_f%ts_pathgeoms, 2)
  
  
  fortran_type_f%ts_lostrans = 0_fpk
  ts_lostrans_ptr => fortran_type_f%ts_lostrans
  transfer_struct_c%ts_lostrans = c_loc(ts_lostrans_ptr(&
    lbound(fortran_type_f%ts_lostrans,1),&
    lbound(fortran_type_f%ts_lostrans,2)))
  inquire(iolength=transfer_struct_c%ts_lostrans_f_byte_size) fortran_type_f%ts_lostrans(&
    lbound(fortran_type_f%ts_lostrans,1),&
    lbound(fortran_type_f%ts_lostrans,2))
#ifdef ifort
  transfer_struct_c%ts_lostrans_f_byte_size = transfer_struct_c%ts_lostrans_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lostrans_f_shapes(1) = size(fortran_type_f%ts_lostrans, 1)
  transfer_struct_c%ts_lostrans_f_shapes(2) = size(fortran_type_f%ts_lostrans, 2)
  
  
  fortran_type_f%ts_layer_mssts = 0_fpk
  ts_layer_mssts_ptr => fortran_type_f%ts_layer_mssts
  transfer_struct_c%ts_layer_mssts = c_loc(ts_layer_mssts_ptr(&
    lbound(fortran_type_f%ts_layer_mssts,1),&
    lbound(fortran_type_f%ts_layer_mssts,2)))
  inquire(iolength=transfer_struct_c%ts_layer_mssts_f_byte_size) fortran_type_f%ts_layer_mssts(&
    lbound(fortran_type_f%ts_layer_mssts,1),&
    lbound(fortran_type_f%ts_layer_mssts,2))
#ifdef ifort
  transfer_struct_c%ts_layer_mssts_f_byte_size = transfer_struct_c%ts_layer_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_layer_mssts_f_shapes(1) = size(fortran_type_f%ts_layer_mssts, 1)
  transfer_struct_c%ts_layer_mssts_f_shapes(2) = size(fortran_type_f%ts_layer_mssts, 2)
  
  
  fortran_type_f%ts_surf_mssts = 0_fpk
  ts_surf_mssts_ptr => fortran_type_f%ts_surf_mssts
  transfer_struct_c%ts_surf_mssts = c_loc(ts_surf_mssts_ptr(&
    lbound(fortran_type_f%ts_surf_mssts,1)))
  inquire(iolength=transfer_struct_c%ts_surf_mssts_f_byte_size) fortran_type_f%ts_surf_mssts(&
    lbound(fortran_type_f%ts_surf_mssts,1))
#ifdef ifort
  transfer_struct_c%ts_surf_mssts_f_byte_size = transfer_struct_c%ts_surf_mssts_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_surf_mssts_f_shapes(1) = size(fortran_type_f%ts_surf_mssts, 1)
  
  
  fortran_type_f%ts_contribs = 0_fpk
  ts_contribs_ptr => fortran_type_f%ts_contribs
  transfer_struct_c%ts_contribs = c_loc(ts_contribs_ptr(&
    lbound(fortran_type_f%ts_contribs,1),&
    lbound(fortran_type_f%ts_contribs,2)))
  inquire(iolength=transfer_struct_c%ts_contribs_f_byte_size) fortran_type_f%ts_contribs(&
    lbound(fortran_type_f%ts_contribs,1),&
    lbound(fortran_type_f%ts_contribs,2))
#ifdef ifort
  transfer_struct_c%ts_contribs_f_byte_size = transfer_struct_c%ts_contribs_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_contribs_f_shapes(1) = size(fortran_type_f%ts_contribs, 1)
  transfer_struct_c%ts_contribs_f_shapes(2) = size(fortran_type_f%ts_contribs, 2)
  
  
  fortran_type_f%ts_fourier_saved = 0
  ts_fourier_saved_ptr => fortran_type_f%ts_fourier_saved
  transfer_struct_c%ts_fourier_saved = c_loc(ts_fourier_saved_ptr(&
    lbound(fortran_type_f%ts_fourier_saved,1)))
  inquire(iolength=transfer_struct_c%ts_fourier_saved_f_byte_size) fortran_type_f%ts_fourier_saved(&
    lbound(fortran_type_f%ts_fourier_saved,1))
#ifdef ifort
  transfer_struct_c%ts_fourier_saved_f_byte_size = transfer_struct_c%ts_fourier_saved_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_fourier_saved_f_shapes(1) = size(fortran_type_f%ts_fourier_saved, 1)
  
  
  fortran_type_f%ts_n_geometries = 0
  ts_n_geometries_ptr => fortran_type_f%ts_n_geometries
  transfer_struct_c%ts_n_geometries = c_loc(ts_n_geometries_ptr)
  inquire(iolength=transfer_struct_c%ts_n_geometries_f_byte_size) fortran_type_f%ts_n_geometries
#ifdef ifort
  transfer_struct_c%ts_n_geometries_f_byte_size = transfer_struct_c%ts_n_geometries_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_solarbeam_boatrans = 0_fpk
  ts_solarbeam_boatrans_ptr => fortran_type_f%ts_solarbeam_boatrans
  transfer_struct_c%ts_solarbeam_boatrans = c_loc(ts_solarbeam_boatrans_ptr(&
    lbound(fortran_type_f%ts_solarbeam_boatrans,1)))
  inquire(iolength=transfer_struct_c%ts_solarbeam_boatrans_f_byte_size) fortran_type_f%ts_solarbeam_boatrans(&
    lbound(fortran_type_f%ts_solarbeam_boatrans,1))
#ifdef ifort
  transfer_struct_c%ts_solarbeam_boatrans_f_byte_size = transfer_struct_c%ts_solarbeam_boatrans_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_solarbeam_boatrans_f_shapes(1) = size(fortran_type_f%ts_solarbeam_boatrans, 1)
  
  
  fortran_type_f%ts_spheralb = 0_fpk
  ts_spheralb_ptr => fortran_type_f%ts_spheralb
  transfer_struct_c%ts_spheralb = c_loc(ts_spheralb_ptr)
  inquire(iolength=transfer_struct_c%ts_spheralb_f_byte_size) fortran_type_f%ts_spheralb
#ifdef ifort
  transfer_struct_c%ts_spheralb_f_byte_size = transfer_struct_c%ts_spheralb_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_trans1_user = 0_fpk
  ts_trans1_user_ptr => fortran_type_f%ts_trans1_user
  transfer_struct_c%ts_trans1_user = c_loc(ts_trans1_user_ptr(&
    lbound(fortran_type_f%ts_trans1_user,1)))
  inquire(iolength=transfer_struct_c%ts_trans1_user_f_byte_size) fortran_type_f%ts_trans1_user(&
    lbound(fortran_type_f%ts_trans1_user,1))
#ifdef ifort
  transfer_struct_c%ts_trans1_user_f_byte_size = transfer_struct_c%ts_trans1_user_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trans1_user_f_shapes(1) = size(fortran_type_f%ts_trans1_user, 1)
  
  
  fortran_type_f%ts_trans1_beam = 0_fpk
  ts_trans1_beam_ptr => fortran_type_f%ts_trans1_beam
  transfer_struct_c%ts_trans1_beam = c_loc(ts_trans1_beam_ptr(&
    lbound(fortran_type_f%ts_trans1_beam,1)))
  inquire(iolength=transfer_struct_c%ts_trans1_beam_f_byte_size) fortran_type_f%ts_trans1_beam(&
    lbound(fortran_type_f%ts_trans1_beam,1))
#ifdef ifort
  transfer_struct_c%ts_trans1_beam_f_byte_size = transfer_struct_c%ts_trans1_beam_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_trans1_beam_f_shapes(1) = size(fortran_type_f%ts_trans1_beam, 1)
  
  
end subroutine lidort_main_outputs_c_init_only

subroutine lidort_main_outputs_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_main_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_main_outputs_c_destroy

subroutine lidort_main_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def_m, only : lidort_main_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_main_outputs), pointer :: fortran_type_f_from
  type(lidort_main_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_intensity = fortran_type_f_from%ts_intensity
  fortran_type_f_to%ts_meani_diffuse = fortran_type_f_from%ts_meani_diffuse
  fortran_type_f_to%ts_flux_diffuse = fortran_type_f_from%ts_flux_diffuse
  fortran_type_f_to%ts_dnmeani_direct = fortran_type_f_from%ts_dnmeani_direct
  fortran_type_f_to%ts_dnflux_direct = fortran_type_f_from%ts_dnflux_direct
  fortran_type_f_to%ts_albmed_user = fortran_type_f_from%ts_albmed_user
  fortran_type_f_to%ts_trnmed_user = fortran_type_f_from%ts_trnmed_user
  fortran_type_f_to%ts_albmed_fluxes = fortran_type_f_from%ts_albmed_fluxes
  fortran_type_f_to%ts_trnmed_fluxes = fortran_type_f_from%ts_trnmed_fluxes
  fortran_type_f_to%ts_planetary_transterm = fortran_type_f_from%ts_planetary_transterm
  fortran_type_f_to%ts_planetary_sbterm = fortran_type_f_from%ts_planetary_sbterm
  fortran_type_f_to%ts_pathgeoms = fortran_type_f_from%ts_pathgeoms
  fortran_type_f_to%ts_lostrans = fortran_type_f_from%ts_lostrans
  fortran_type_f_to%ts_layer_mssts = fortran_type_f_from%ts_layer_mssts
  fortran_type_f_to%ts_surf_mssts = fortran_type_f_from%ts_surf_mssts
  fortran_type_f_to%ts_contribs = fortran_type_f_from%ts_contribs
  fortran_type_f_to%ts_fourier_saved = fortran_type_f_from%ts_fourier_saved
  fortran_type_f_to%ts_n_geometries = fortran_type_f_from%ts_n_geometries
  fortran_type_f_to%ts_solarbeam_boatrans = fortran_type_f_from%ts_solarbeam_boatrans
  fortran_type_f_to%ts_spheralb = fortran_type_f_from%ts_spheralb
  fortran_type_f_to%ts_trans1_user = fortran_type_f_from%ts_trans1_user
  fortran_type_f_to%ts_trans1_beam = fortran_type_f_from%ts_trans1_beam
  

end subroutine lidort_main_outputs_c_copy

! Links to type: "lidort_wladjusted_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_wladjusted_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_wladjusted_outputs

  type(lidort_wladjusted_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_wladjusted_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_wladjusted_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_wladjusted_outputs_c_alloc_init

! Links to type: "lidort_wladjusted_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_wladjusted_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m
  use lidort_pars_m

  type(lidort_wladjusted_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_wladjusted_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:), pointer :: ts_wladjusted_isotropic_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_wladjusted_direct_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_wladjusted_f_ords_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_wladjusted_f_user_0_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_wladjusted_isotropic = 0_fpk
  ts_wladjusted_isotropic_ptr => fortran_type_f%ts_wladjusted_isotropic
  transfer_struct_c%ts_wladjusted_isotropic = c_loc(ts_wladjusted_isotropic_ptr(&
    lbound(fortran_type_f%ts_wladjusted_isotropic,1)))
  inquire(iolength=transfer_struct_c%ts_wladjusted_isotropic_f_byte_size) fortran_type_f%ts_wladjusted_isotropic(&
    lbound(fortran_type_f%ts_wladjusted_isotropic,1))
#ifdef ifort
  transfer_struct_c%ts_wladjusted_isotropic_f_byte_size = transfer_struct_c%ts_wladjusted_isotropic_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_wladjusted_isotropic_f_shapes(1) = size(fortran_type_f%ts_wladjusted_isotropic, 1)
  
  
  fortran_type_f%ts_wladjusted_direct = 0_fpk
  ts_wladjusted_direct_ptr => fortran_type_f%ts_wladjusted_direct
  transfer_struct_c%ts_wladjusted_direct = c_loc(ts_wladjusted_direct_ptr(&
    lbound(fortran_type_f%ts_wladjusted_direct,1),&
    lbound(fortran_type_f%ts_wladjusted_direct,2),&
    lbound(fortran_type_f%ts_wladjusted_direct,3)))
  inquire(iolength=transfer_struct_c%ts_wladjusted_direct_f_byte_size) fortran_type_f%ts_wladjusted_direct(&
    lbound(fortran_type_f%ts_wladjusted_direct,1),&
    lbound(fortran_type_f%ts_wladjusted_direct,2),&
    lbound(fortran_type_f%ts_wladjusted_direct,3))
#ifdef ifort
  transfer_struct_c%ts_wladjusted_direct_f_byte_size = transfer_struct_c%ts_wladjusted_direct_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_wladjusted_direct_f_shapes(1) = size(fortran_type_f%ts_wladjusted_direct, 1)
  transfer_struct_c%ts_wladjusted_direct_f_shapes(2) = size(fortran_type_f%ts_wladjusted_direct, 2)
  transfer_struct_c%ts_wladjusted_direct_f_shapes(3) = size(fortran_type_f%ts_wladjusted_direct, 3)
  
  
  fortran_type_f%ts_wladjusted_f_ords_0 = 0_fpk
  ts_wladjusted_f_ords_0_ptr => fortran_type_f%ts_wladjusted_f_ords_0
  transfer_struct_c%ts_wladjusted_f_ords_0 = c_loc(ts_wladjusted_f_ords_0_ptr(&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,1),&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,2),&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,3)))
  inquire(iolength=transfer_struct_c%ts_wladjusted_f_ords_0_f_byte_size) fortran_type_f%ts_wladjusted_f_ords_0(&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,1),&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,2),&
    lbound(fortran_type_f%ts_wladjusted_f_ords_0,3))
#ifdef ifort
  transfer_struct_c%ts_wladjusted_f_ords_0_f_byte_size = transfer_struct_c%ts_wladjusted_f_ords_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_wladjusted_f_ords_0_f_shapes(1) = size(fortran_type_f%ts_wladjusted_f_ords_0, 1)
  transfer_struct_c%ts_wladjusted_f_ords_0_f_shapes(2) = size(fortran_type_f%ts_wladjusted_f_ords_0, 2)
  transfer_struct_c%ts_wladjusted_f_ords_0_f_shapes(3) = size(fortran_type_f%ts_wladjusted_f_ords_0, 3)
  
  
  fortran_type_f%ts_wladjusted_f_user_0 = 0_fpk
  ts_wladjusted_f_user_0_ptr => fortran_type_f%ts_wladjusted_f_user_0
  transfer_struct_c%ts_wladjusted_f_user_0 = c_loc(ts_wladjusted_f_user_0_ptr(&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,1),&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,2),&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,3)))
  inquire(iolength=transfer_struct_c%ts_wladjusted_f_user_0_f_byte_size) fortran_type_f%ts_wladjusted_f_user_0(&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,1),&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,2),&
    lbound(fortran_type_f%ts_wladjusted_f_user_0,3))
#ifdef ifort
  transfer_struct_c%ts_wladjusted_f_user_0_f_byte_size = transfer_struct_c%ts_wladjusted_f_user_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_wladjusted_f_user_0_f_shapes(1) = size(fortran_type_f%ts_wladjusted_f_user_0, 1)
  transfer_struct_c%ts_wladjusted_f_user_0_f_shapes(2) = size(fortran_type_f%ts_wladjusted_f_user_0, 2)
  transfer_struct_c%ts_wladjusted_f_user_0_f_shapes(3) = size(fortran_type_f%ts_wladjusted_f_user_0, 3)
  
  
end subroutine lidort_wladjusted_outputs_c_init_only

subroutine lidort_wladjusted_outputs_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_wladjusted_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_wladjusted_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_wladjusted_outputs_c_destroy

subroutine lidort_wladjusted_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def_m, only : lidort_wladjusted_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_wladjusted_outputs), pointer :: fortran_type_f_from
  type(lidort_wladjusted_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_wladjusted_isotropic = fortran_type_f_from%ts_wladjusted_isotropic
  fortran_type_f_to%ts_wladjusted_direct = fortran_type_f_from%ts_wladjusted_direct
  fortran_type_f_to%ts_wladjusted_f_ords_0 = fortran_type_f_from%ts_wladjusted_f_ords_0
  fortran_type_f_to%ts_wladjusted_f_user_0 = fortran_type_f_from%ts_wladjusted_f_user_0
  

end subroutine lidort_wladjusted_outputs_c_copy

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(lidort_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_exception_handling_c_alloc_init

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m
  use lidort_pars_m

  type(lidort_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_status_inputcheck_ptr
  integer(c_int), pointer :: ts_ncheckmessages_ptr
  integer(c_int), pointer :: ts_status_calculation_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_status_inputcheck = 0
  ts_status_inputcheck_ptr => fortran_type_f%ts_status_inputcheck
  transfer_struct_c%ts_status_inputcheck = c_loc(ts_status_inputcheck_ptr)
  inquire(iolength=transfer_struct_c%ts_status_inputcheck_f_byte_size) fortran_type_f%ts_status_inputcheck
#ifdef ifort
  transfer_struct_c%ts_status_inputcheck_f_byte_size = transfer_struct_c%ts_status_inputcheck_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_ncheckmessages = 0
  ts_ncheckmessages_ptr => fortran_type_f%ts_ncheckmessages
  transfer_struct_c%ts_ncheckmessages = c_loc(ts_ncheckmessages_ptr)
  inquire(iolength=transfer_struct_c%ts_ncheckmessages_f_byte_size) fortran_type_f%ts_ncheckmessages
#ifdef ifort
  transfer_struct_c%ts_ncheckmessages_f_byte_size = transfer_struct_c%ts_ncheckmessages_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_checkmessages = ''
  transfer_struct_c%ts_checkmessages_f_len = len(fortran_type_f%ts_checkmessages)
  transfer_struct_c%ts_checkmessages_f_shapes(1) = size(fortran_type_f%ts_checkmessages, 1)
  
  fortran_type_f%ts_actions = ''
  transfer_struct_c%ts_actions_f_len = len(fortran_type_f%ts_actions)
  transfer_struct_c%ts_actions_f_shapes(1) = size(fortran_type_f%ts_actions, 1)
  
  
  fortran_type_f%ts_status_calculation = 0
  ts_status_calculation_ptr => fortran_type_f%ts_status_calculation
  transfer_struct_c%ts_status_calculation = c_loc(ts_status_calculation_ptr)
  inquire(iolength=transfer_struct_c%ts_status_calculation_f_byte_size) fortran_type_f%ts_status_calculation
#ifdef ifort
  transfer_struct_c%ts_status_calculation_f_byte_size = transfer_struct_c%ts_status_calculation_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_message = ''
  transfer_struct_c%ts_message_f_len = len(fortran_type_f%ts_message)
  
  fortran_type_f%ts_trace_1 = ''
  transfer_struct_c%ts_trace_1_f_len = len(fortran_type_f%ts_trace_1)
  
  fortran_type_f%ts_trace_2 = ''
  transfer_struct_c%ts_trace_2_f_len = len(fortran_type_f%ts_trace_2)
  
  fortran_type_f%ts_trace_3 = ''
  transfer_struct_c%ts_trace_3_f_len = len(fortran_type_f%ts_trace_3)
  
  
end subroutine lidort_exception_handling_c_init_only

subroutine lidort_exception_handling_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_exception_handling_c_destroy

subroutine lidort_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_exception_handling), pointer :: fortran_type_f_from
  type(lidort_exception_handling), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_status_inputcheck = fortran_type_f_from%ts_status_inputcheck
  fortran_type_f_to%ts_ncheckmessages = fortran_type_f_from%ts_ncheckmessages
  fortran_type_f_to%ts_checkmessages = fortran_type_f_from%ts_checkmessages
  fortran_type_f_to%ts_actions = fortran_type_f_from%ts_actions
  fortran_type_f_to%ts_status_calculation = fortran_type_f_from%ts_status_calculation
  fortran_type_f_to%ts_message = fortran_type_f_from%ts_message
  fortran_type_f_to%ts_trace_1 = fortran_type_f_from%ts_trace_1
  fortran_type_f_to%ts_trace_2 = fortran_type_f_from%ts_trace_2
  fortran_type_f_to%ts_trace_3 = fortran_type_f_from%ts_trace_3
  

end subroutine lidort_exception_handling_c_copy

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_input_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_input_exception_handling

  type(lidort_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_input_exception_handling_c_alloc_init

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m
  use lidort_pars_m

  type(lidort_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_status_inputread_ptr
  integer(c_int), pointer :: ts_ninputmessages_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_status_inputread = 0
  ts_status_inputread_ptr => fortran_type_f%ts_status_inputread
  transfer_struct_c%ts_status_inputread = c_loc(ts_status_inputread_ptr)
  inquire(iolength=transfer_struct_c%ts_status_inputread_f_byte_size) fortran_type_f%ts_status_inputread
#ifdef ifort
  transfer_struct_c%ts_status_inputread_f_byte_size = transfer_struct_c%ts_status_inputread_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_ninputmessages = 0
  ts_ninputmessages_ptr => fortran_type_f%ts_ninputmessages
  transfer_struct_c%ts_ninputmessages = c_loc(ts_ninputmessages_ptr)
  inquire(iolength=transfer_struct_c%ts_ninputmessages_f_byte_size) fortran_type_f%ts_ninputmessages
#ifdef ifort
  transfer_struct_c%ts_ninputmessages_f_byte_size = transfer_struct_c%ts_ninputmessages_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_inputmessages = ''
  transfer_struct_c%ts_inputmessages_f_len = len(fortran_type_f%ts_inputmessages)
  transfer_struct_c%ts_inputmessages_f_shapes(1) = size(fortran_type_f%ts_inputmessages, 1)
  
  fortran_type_f%ts_inputactions = ''
  transfer_struct_c%ts_inputactions_f_len = len(fortran_type_f%ts_inputactions)
  transfer_struct_c%ts_inputactions_f_shapes(1) = size(fortran_type_f%ts_inputactions, 1)
  
  
end subroutine lidort_input_exception_handling_c_init_only

subroutine lidort_input_exception_handling_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_input_exception_handling_c_destroy

subroutine lidort_input_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def_m, only : lidort_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_input_exception_handling), pointer :: fortran_type_f_from
  type(lidort_input_exception_handling), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_status_inputread = fortran_type_f_from%ts_status_inputread
  fortran_type_f_to%ts_ninputmessages = fortran_type_f_from%ts_ninputmessages
  fortran_type_f_to%ts_inputmessages = fortran_type_f_from%ts_inputmessages
  fortran_type_f_to%ts_inputactions = fortran_type_f_from%ts_inputactions
  

end subroutine lidort_input_exception_handling_c_copy

! Links to type: "lidort_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_outputs

  type(lidort_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_outputs_c_alloc_init

! Links to type: "lidort_outputs" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def_m
  use lidort_pars_m

  type(lidort_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  type(lidort_main_outputs), pointer :: main_ptr
  type(lidort_wladjusted_outputs), pointer :: wlout_ptr
  type(lidort_exception_handling), pointer :: status_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  main_ptr => fortran_type_f%main
  transfer_struct_c%main = c_loc(main_ptr)
  inquire(iolength=transfer_struct_c%main_f_byte_size) fortran_type_f%main
#ifdef ifort
  transfer_struct_c%main_f_byte_size = transfer_struct_c%main_f_byte_size * 4
#endif
  
  
  wlout_ptr => fortran_type_f%wlout
  transfer_struct_c%wlout = c_loc(wlout_ptr)
  inquire(iolength=transfer_struct_c%wlout_f_byte_size) fortran_type_f%wlout
#ifdef ifort
  transfer_struct_c%wlout_f_byte_size = transfer_struct_c%wlout_f_byte_size * 4
#endif
  
  
  status_ptr => fortran_type_f%status
  transfer_struct_c%status = c_loc(status_ptr)
  inquire(iolength=transfer_struct_c%status_f_byte_size) fortran_type_f%status
#ifdef ifort
  transfer_struct_c%status_f_byte_size = transfer_struct_c%status_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_outputs_c_init_only

subroutine lidort_outputs_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def_m, only : lidort_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_outputs_c_destroy

subroutine lidort_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def_m, only : lidort_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_outputs), pointer :: fortran_type_f_from
  type(lidort_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%main = fortran_type_f_from%main
  fortran_type_f_to%wlout = fortran_type_f_from%wlout
  fortran_type_f_to%status = fortran_type_f_from%status
  

end subroutine lidort_outputs_c_copy

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def_m" in file: "lidort_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_sup_brdf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_brdf_def_m, only : lidort_sup_brdf

  type(lidort_sup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_brdf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_brdf_c_alloc_init

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def_m" in file: "lidort_sup_brdf_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_brdf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_brdf_def_m
  use lidort_pars_m

  type(lidort_sup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_user_brdf_f_ptr
  real(c_double), dimension(:), pointer :: ts_emissivity_ptr
  real(c_double), dimension(:), pointer :: ts_user_emissivity_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_exactdb_brdfunc = 0_fpk
  ts_exactdb_brdfunc_ptr => fortran_type_f%ts_exactdb_brdfunc
  transfer_struct_c%ts_exactdb_brdfunc = c_loc(ts_exactdb_brdfunc_ptr(&
    lbound(fortran_type_f%ts_exactdb_brdfunc,1),&
    lbound(fortran_type_f%ts_exactdb_brdfunc,2),&
    lbound(fortran_type_f%ts_exactdb_brdfunc,3)))
  inquire(iolength=transfer_struct_c%ts_exactdb_brdfunc_f_byte_size) fortran_type_f%ts_exactdb_brdfunc(&
    lbound(fortran_type_f%ts_exactdb_brdfunc,1),&
    lbound(fortran_type_f%ts_exactdb_brdfunc,2),&
    lbound(fortran_type_f%ts_exactdb_brdfunc,3))
#ifdef ifort
  transfer_struct_c%ts_exactdb_brdfunc_f_byte_size = transfer_struct_c%ts_exactdb_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_exactdb_brdfunc_f_shapes(1) = size(fortran_type_f%ts_exactdb_brdfunc, 1)
  transfer_struct_c%ts_exactdb_brdfunc_f_shapes(2) = size(fortran_type_f%ts_exactdb_brdfunc, 2)
  transfer_struct_c%ts_exactdb_brdfunc_f_shapes(3) = size(fortran_type_f%ts_exactdb_brdfunc, 3)
  
  
  fortran_type_f%ts_brdf_f_0 = 0_fpk
  ts_brdf_f_0_ptr => fortran_type_f%ts_brdf_f_0
  transfer_struct_c%ts_brdf_f_0 = c_loc(ts_brdf_f_0_ptr(&
    lbound(fortran_type_f%ts_brdf_f_0,1),&
    lbound(fortran_type_f%ts_brdf_f_0,2),&
    lbound(fortran_type_f%ts_brdf_f_0,3)))
  inquire(iolength=transfer_struct_c%ts_brdf_f_0_f_byte_size) fortran_type_f%ts_brdf_f_0(&
    lbound(fortran_type_f%ts_brdf_f_0,1),&
    lbound(fortran_type_f%ts_brdf_f_0,2),&
    lbound(fortran_type_f%ts_brdf_f_0,3))
#ifdef ifort
  transfer_struct_c%ts_brdf_f_0_f_byte_size = transfer_struct_c%ts_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_brdf_f_0_f_shapes(1) = size(fortran_type_f%ts_brdf_f_0, 1)
  transfer_struct_c%ts_brdf_f_0_f_shapes(2) = size(fortran_type_f%ts_brdf_f_0, 2)
  transfer_struct_c%ts_brdf_f_0_f_shapes(3) = size(fortran_type_f%ts_brdf_f_0, 3)
  
  
  fortran_type_f%ts_brdf_f = 0_fpk
  ts_brdf_f_ptr => fortran_type_f%ts_brdf_f
  transfer_struct_c%ts_brdf_f = c_loc(ts_brdf_f_ptr(&
    lbound(fortran_type_f%ts_brdf_f,1),&
    lbound(fortran_type_f%ts_brdf_f,2),&
    lbound(fortran_type_f%ts_brdf_f,3)))
  inquire(iolength=transfer_struct_c%ts_brdf_f_f_byte_size) fortran_type_f%ts_brdf_f(&
    lbound(fortran_type_f%ts_brdf_f,1),&
    lbound(fortran_type_f%ts_brdf_f,2),&
    lbound(fortran_type_f%ts_brdf_f,3))
#ifdef ifort
  transfer_struct_c%ts_brdf_f_f_byte_size = transfer_struct_c%ts_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_brdf_f_f_shapes(1) = size(fortran_type_f%ts_brdf_f, 1)
  transfer_struct_c%ts_brdf_f_f_shapes(2) = size(fortran_type_f%ts_brdf_f, 2)
  transfer_struct_c%ts_brdf_f_f_shapes(3) = size(fortran_type_f%ts_brdf_f, 3)
  
  
  fortran_type_f%ts_user_brdf_f_0 = 0_fpk
  ts_user_brdf_f_0_ptr => fortran_type_f%ts_user_brdf_f_0
  transfer_struct_c%ts_user_brdf_f_0 = c_loc(ts_user_brdf_f_0_ptr(&
    lbound(fortran_type_f%ts_user_brdf_f_0,1),&
    lbound(fortran_type_f%ts_user_brdf_f_0,2),&
    lbound(fortran_type_f%ts_user_brdf_f_0,3)))
  inquire(iolength=transfer_struct_c%ts_user_brdf_f_0_f_byte_size) fortran_type_f%ts_user_brdf_f_0(&
    lbound(fortran_type_f%ts_user_brdf_f_0,1),&
    lbound(fortran_type_f%ts_user_brdf_f_0,2),&
    lbound(fortran_type_f%ts_user_brdf_f_0,3))
#ifdef ifort
  transfer_struct_c%ts_user_brdf_f_0_f_byte_size = transfer_struct_c%ts_user_brdf_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_brdf_f_0_f_shapes(1) = size(fortran_type_f%ts_user_brdf_f_0, 1)
  transfer_struct_c%ts_user_brdf_f_0_f_shapes(2) = size(fortran_type_f%ts_user_brdf_f_0, 2)
  transfer_struct_c%ts_user_brdf_f_0_f_shapes(3) = size(fortran_type_f%ts_user_brdf_f_0, 3)
  
  
  fortran_type_f%ts_user_brdf_f = 0_fpk
  ts_user_brdf_f_ptr => fortran_type_f%ts_user_brdf_f
  transfer_struct_c%ts_user_brdf_f = c_loc(ts_user_brdf_f_ptr(&
    lbound(fortran_type_f%ts_user_brdf_f,1),&
    lbound(fortran_type_f%ts_user_brdf_f,2),&
    lbound(fortran_type_f%ts_user_brdf_f,3)))
  inquire(iolength=transfer_struct_c%ts_user_brdf_f_f_byte_size) fortran_type_f%ts_user_brdf_f(&
    lbound(fortran_type_f%ts_user_brdf_f,1),&
    lbound(fortran_type_f%ts_user_brdf_f,2),&
    lbound(fortran_type_f%ts_user_brdf_f,3))
#ifdef ifort
  transfer_struct_c%ts_user_brdf_f_f_byte_size = transfer_struct_c%ts_user_brdf_f_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_brdf_f_f_shapes(1) = size(fortran_type_f%ts_user_brdf_f, 1)
  transfer_struct_c%ts_user_brdf_f_f_shapes(2) = size(fortran_type_f%ts_user_brdf_f, 2)
  transfer_struct_c%ts_user_brdf_f_f_shapes(3) = size(fortran_type_f%ts_user_brdf_f, 3)
  
  
  fortran_type_f%ts_emissivity = 0_fpk
  ts_emissivity_ptr => fortran_type_f%ts_emissivity
  transfer_struct_c%ts_emissivity = c_loc(ts_emissivity_ptr(&
    lbound(fortran_type_f%ts_emissivity,1)))
  inquire(iolength=transfer_struct_c%ts_emissivity_f_byte_size) fortran_type_f%ts_emissivity(&
    lbound(fortran_type_f%ts_emissivity,1))
#ifdef ifort
  transfer_struct_c%ts_emissivity_f_byte_size = transfer_struct_c%ts_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_emissivity_f_shapes(1) = size(fortran_type_f%ts_emissivity, 1)
  
  
  fortran_type_f%ts_user_emissivity = 0_fpk
  ts_user_emissivity_ptr => fortran_type_f%ts_user_emissivity
  transfer_struct_c%ts_user_emissivity = c_loc(ts_user_emissivity_ptr(&
    lbound(fortran_type_f%ts_user_emissivity,1)))
  inquire(iolength=transfer_struct_c%ts_user_emissivity_f_byte_size) fortran_type_f%ts_user_emissivity(&
    lbound(fortran_type_f%ts_user_emissivity,1))
#ifdef ifort
  transfer_struct_c%ts_user_emissivity_f_byte_size = transfer_struct_c%ts_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_emissivity_f_shapes(1) = size(fortran_type_f%ts_user_emissivity, 1)
  
  
end subroutine lidort_sup_brdf_c_init_only

subroutine lidort_sup_brdf_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_brdf_def_m, only : lidort_sup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_brdf_c_destroy

subroutine lidort_sup_brdf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_brdf_def_m, only : lidort_sup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_sup_brdf), pointer :: fortran_type_f_from
  type(lidort_sup_brdf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_exactdb_brdfunc = fortran_type_f_from%ts_exactdb_brdfunc
  fortran_type_f_to%ts_brdf_f_0 = fortran_type_f_from%ts_brdf_f_0
  fortran_type_f_to%ts_brdf_f = fortran_type_f_from%ts_brdf_f
  fortran_type_f_to%ts_user_brdf_f_0 = fortran_type_f_from%ts_user_brdf_f_0
  fortran_type_f_to%ts_user_brdf_f = fortran_type_f_from%ts_user_brdf_f
  fortran_type_f_to%ts_emissivity = fortran_type_f_from%ts_emissivity
  fortran_type_f_to%ts_user_emissivity = fortran_type_f_from%ts_user_emissivity
  

end subroutine lidort_sup_brdf_c_copy

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def_m" in file: "lidort_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_sup_sleave_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_sleave_def_m, only : lidort_sup_sleave

  type(lidort_sup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_sleave_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_sleave_c_alloc_init

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def_m" in file: "lidort_sup_sleave_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_sleave_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_sleave_def_m
  use lidort_pars_m

  type(lidort_sup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  real(c_double), dimension(:), pointer :: ts_slterm_isotropic_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_slterm_userangles_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_slterm_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_user_slterm_f_0_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_slterm_isotropic = 0_fpk
  ts_slterm_isotropic_ptr => fortran_type_f%ts_slterm_isotropic
  transfer_struct_c%ts_slterm_isotropic = c_loc(ts_slterm_isotropic_ptr(&
    lbound(fortran_type_f%ts_slterm_isotropic,1)))
  inquire(iolength=transfer_struct_c%ts_slterm_isotropic_f_byte_size) fortran_type_f%ts_slterm_isotropic(&
    lbound(fortran_type_f%ts_slterm_isotropic,1))
#ifdef ifort
  transfer_struct_c%ts_slterm_isotropic_f_byte_size = transfer_struct_c%ts_slterm_isotropic_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_slterm_isotropic_f_shapes(1) = size(fortran_type_f%ts_slterm_isotropic, 1)
  
  
  fortran_type_f%ts_slterm_userangles = 0_fpk
  ts_slterm_userangles_ptr => fortran_type_f%ts_slterm_userangles
  transfer_struct_c%ts_slterm_userangles = c_loc(ts_slterm_userangles_ptr(&
    lbound(fortran_type_f%ts_slterm_userangles,1),&
    lbound(fortran_type_f%ts_slterm_userangles,2),&
    lbound(fortran_type_f%ts_slterm_userangles,3)))
  inquire(iolength=transfer_struct_c%ts_slterm_userangles_f_byte_size) fortran_type_f%ts_slterm_userangles(&
    lbound(fortran_type_f%ts_slterm_userangles,1),&
    lbound(fortran_type_f%ts_slterm_userangles,2),&
    lbound(fortran_type_f%ts_slterm_userangles,3))
#ifdef ifort
  transfer_struct_c%ts_slterm_userangles_f_byte_size = transfer_struct_c%ts_slterm_userangles_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_slterm_userangles_f_shapes(1) = size(fortran_type_f%ts_slterm_userangles, 1)
  transfer_struct_c%ts_slterm_userangles_f_shapes(2) = size(fortran_type_f%ts_slterm_userangles, 2)
  transfer_struct_c%ts_slterm_userangles_f_shapes(3) = size(fortran_type_f%ts_slterm_userangles, 3)
  
  
  fortran_type_f%ts_slterm_f_0 = 0_fpk
  ts_slterm_f_0_ptr => fortran_type_f%ts_slterm_f_0
  transfer_struct_c%ts_slterm_f_0 = c_loc(ts_slterm_f_0_ptr(&
    lbound(fortran_type_f%ts_slterm_f_0,1),&
    lbound(fortran_type_f%ts_slterm_f_0,2),&
    lbound(fortran_type_f%ts_slterm_f_0,3)))
  inquire(iolength=transfer_struct_c%ts_slterm_f_0_f_byte_size) fortran_type_f%ts_slterm_f_0(&
    lbound(fortran_type_f%ts_slterm_f_0,1),&
    lbound(fortran_type_f%ts_slterm_f_0,2),&
    lbound(fortran_type_f%ts_slterm_f_0,3))
#ifdef ifort
  transfer_struct_c%ts_slterm_f_0_f_byte_size = transfer_struct_c%ts_slterm_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_slterm_f_0_f_shapes(1) = size(fortran_type_f%ts_slterm_f_0, 1)
  transfer_struct_c%ts_slterm_f_0_f_shapes(2) = size(fortran_type_f%ts_slterm_f_0, 2)
  transfer_struct_c%ts_slterm_f_0_f_shapes(3) = size(fortran_type_f%ts_slterm_f_0, 3)
  
  
  fortran_type_f%ts_user_slterm_f_0 = 0_fpk
  ts_user_slterm_f_0_ptr => fortran_type_f%ts_user_slterm_f_0
  transfer_struct_c%ts_user_slterm_f_0 = c_loc(ts_user_slterm_f_0_ptr(&
    lbound(fortran_type_f%ts_user_slterm_f_0,1),&
    lbound(fortran_type_f%ts_user_slterm_f_0,2),&
    lbound(fortran_type_f%ts_user_slterm_f_0,3)))
  inquire(iolength=transfer_struct_c%ts_user_slterm_f_0_f_byte_size) fortran_type_f%ts_user_slterm_f_0(&
    lbound(fortran_type_f%ts_user_slterm_f_0,1),&
    lbound(fortran_type_f%ts_user_slterm_f_0,2),&
    lbound(fortran_type_f%ts_user_slterm_f_0,3))
#ifdef ifort
  transfer_struct_c%ts_user_slterm_f_0_f_byte_size = transfer_struct_c%ts_user_slterm_f_0_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_slterm_f_0_f_shapes(1) = size(fortran_type_f%ts_user_slterm_f_0, 1)
  transfer_struct_c%ts_user_slterm_f_0_f_shapes(2) = size(fortran_type_f%ts_user_slterm_f_0, 2)
  transfer_struct_c%ts_user_slterm_f_0_f_shapes(3) = size(fortran_type_f%ts_user_slterm_f_0, 3)
  
  
end subroutine lidort_sup_sleave_c_init_only

subroutine lidort_sup_sleave_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_sleave_def_m, only : lidort_sup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_sleave_c_destroy

subroutine lidort_sup_sleave_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_sleave_def_m, only : lidort_sup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_sup_sleave), pointer :: fortran_type_f_from
  type(lidort_sup_sleave), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_slterm_isotropic = fortran_type_f_from%ts_slterm_isotropic
  fortran_type_f_to%ts_slterm_userangles = fortran_type_f_from%ts_slterm_userangles
  fortran_type_f_to%ts_slterm_f_0 = fortran_type_f_from%ts_slterm_f_0
  fortran_type_f_to%ts_user_slterm_f_0 = fortran_type_f_from%ts_user_slterm_f_0
  

end subroutine lidort_sup_sleave_c_copy

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def_m" in file: "lidort_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_sup_ss_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_ss_def_m, only : lidort_sup_ss

  type(lidort_sup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_ss_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_ss_c_alloc_init

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def_m" in file: "lidort_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_ss_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_ss_def_m
  use lidort_pars_m

  type(lidort_sup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_intensity_ss_ptr
  real(c_double), dimension(:,:), pointer :: ts_intensity_db_ptr
  real(c_double), dimension(:,:), pointer :: ts_contribs_ss_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_intensity_ss = 0_fpk
  ts_intensity_ss_ptr => fortran_type_f%ts_intensity_ss
  transfer_struct_c%ts_intensity_ss = c_loc(ts_intensity_ss_ptr(&
    lbound(fortran_type_f%ts_intensity_ss,1),&
    lbound(fortran_type_f%ts_intensity_ss,2),&
    lbound(fortran_type_f%ts_intensity_ss,3)))
  inquire(iolength=transfer_struct_c%ts_intensity_ss_f_byte_size) fortran_type_f%ts_intensity_ss(&
    lbound(fortran_type_f%ts_intensity_ss,1),&
    lbound(fortran_type_f%ts_intensity_ss,2),&
    lbound(fortran_type_f%ts_intensity_ss,3))
#ifdef ifort
  transfer_struct_c%ts_intensity_ss_f_byte_size = transfer_struct_c%ts_intensity_ss_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_intensity_ss_f_shapes(1) = size(fortran_type_f%ts_intensity_ss, 1)
  transfer_struct_c%ts_intensity_ss_f_shapes(2) = size(fortran_type_f%ts_intensity_ss, 2)
  transfer_struct_c%ts_intensity_ss_f_shapes(3) = size(fortran_type_f%ts_intensity_ss, 3)
  
  
  fortran_type_f%ts_intensity_db = 0_fpk
  ts_intensity_db_ptr => fortran_type_f%ts_intensity_db
  transfer_struct_c%ts_intensity_db = c_loc(ts_intensity_db_ptr(&
    lbound(fortran_type_f%ts_intensity_db,1),&
    lbound(fortran_type_f%ts_intensity_db,2)))
  inquire(iolength=transfer_struct_c%ts_intensity_db_f_byte_size) fortran_type_f%ts_intensity_db(&
    lbound(fortran_type_f%ts_intensity_db,1),&
    lbound(fortran_type_f%ts_intensity_db,2))
#ifdef ifort
  transfer_struct_c%ts_intensity_db_f_byte_size = transfer_struct_c%ts_intensity_db_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_intensity_db_f_shapes(1) = size(fortran_type_f%ts_intensity_db, 1)
  transfer_struct_c%ts_intensity_db_f_shapes(2) = size(fortran_type_f%ts_intensity_db, 2)
  
  
  fortran_type_f%ts_contribs_ss = 0_fpk
  ts_contribs_ss_ptr => fortran_type_f%ts_contribs_ss
  transfer_struct_c%ts_contribs_ss = c_loc(ts_contribs_ss_ptr(&
    lbound(fortran_type_f%ts_contribs_ss,1),&
    lbound(fortran_type_f%ts_contribs_ss,2)))
  inquire(iolength=transfer_struct_c%ts_contribs_ss_f_byte_size) fortran_type_f%ts_contribs_ss(&
    lbound(fortran_type_f%ts_contribs_ss,1),&
    lbound(fortran_type_f%ts_contribs_ss,2))
#ifdef ifort
  transfer_struct_c%ts_contribs_ss_f_byte_size = transfer_struct_c%ts_contribs_ss_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_contribs_ss_f_shapes(1) = size(fortran_type_f%ts_contribs_ss, 1)
  transfer_struct_c%ts_contribs_ss_f_shapes(2) = size(fortran_type_f%ts_contribs_ss, 2)
  
  
end subroutine lidort_sup_ss_c_init_only

subroutine lidort_sup_ss_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_ss_def_m, only : lidort_sup_ss

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_ss_c_destroy

subroutine lidort_sup_ss_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_ss_def_m, only : lidort_sup_ss

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_sup_ss), pointer :: fortran_type_f_from
  type(lidort_sup_ss), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_intensity_ss = fortran_type_f_from%ts_intensity_ss
  fortran_type_f_to%ts_intensity_db = fortran_type_f_from%ts_intensity_db
  fortran_type_f_to%ts_contribs_ss = fortran_type_f_from%ts_contribs_ss
  

end subroutine lidort_sup_ss_c_copy

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def_m" in file: "lidort_sup_def.F90"
! Allocs and initializes type
subroutine lidort_sup_inout_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_inout_def_m, only : lidort_sup_inout

  type(lidort_sup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_inout_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_inout_c_alloc_init

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def_m" in file: "lidort_sup_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_inout_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_inout_def_m
  use lidort_sup_brdf_def_m
  use lidort_sup_ss_def_m
  use lidort_sup_sleave_def_m

  type(lidort_sup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  type(lidort_sup_brdf), pointer :: brdf_ptr
  type(lidort_sup_ss), pointer :: ss_ptr
  type(lidort_sup_sleave), pointer :: sleave_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  brdf_ptr => fortran_type_f%brdf
  transfer_struct_c%brdf = c_loc(brdf_ptr)
  inquire(iolength=transfer_struct_c%brdf_f_byte_size) fortran_type_f%brdf
#ifdef ifort
  transfer_struct_c%brdf_f_byte_size = transfer_struct_c%brdf_f_byte_size * 4
#endif
  
  
  ss_ptr => fortran_type_f%ss
  transfer_struct_c%ss = c_loc(ss_ptr)
  inquire(iolength=transfer_struct_c%ss_f_byte_size) fortran_type_f%ss
#ifdef ifort
  transfer_struct_c%ss_f_byte_size = transfer_struct_c%ss_f_byte_size * 4
#endif
  
  
  sleave_ptr => fortran_type_f%sleave
  transfer_struct_c%sleave = c_loc(sleave_ptr)
  inquire(iolength=transfer_struct_c%sleave_f_byte_size) fortran_type_f%sleave
#ifdef ifort
  transfer_struct_c%sleave_f_byte_size = transfer_struct_c%sleave_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_sup_inout_c_init_only

subroutine lidort_sup_inout_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_inout_def_m, only : lidort_sup_inout

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_inout_c_destroy

subroutine lidort_sup_inout_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_inout_def_m, only : lidort_sup_inout

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_sup_inout), pointer :: fortran_type_f_from
  type(lidort_sup_inout), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%brdf = fortran_type_f_from%brdf
  fortran_type_f_to%ss = fortran_type_f_from%ss
  fortran_type_f_to%sleave = fortran_type_f_from%sleave
  

end subroutine lidort_sup_inout_c_copy

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_boolean_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_boolean

  type(lidort_fixed_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_boolean_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_boolean_c_alloc_init

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_boolean_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_fullrad_mode_ptr
  logical(kind=4), pointer :: ts_do_thermal_emission_ptr
  logical(kind=4), pointer :: ts_do_surface_emission_ptr
  logical(kind=4), pointer :: ts_do_plane_parallel_ptr
  logical(kind=4), pointer :: ts_do_brdf_surface_ptr
  logical(kind=4), pointer :: ts_do_upwelling_ptr
  logical(kind=4), pointer :: ts_do_dnwelling_ptr
  logical(kind=4), pointer :: ts_do_toa_contribs_ptr
  logical(kind=4), pointer :: ts_do_surface_leaving_ptr
  logical(kind=4), pointer :: ts_do_sl_isotropic_ptr
  logical(kind=4), pointer :: ts_do_water_leaving_ptr
  logical(kind=4), pointer :: ts_do_fluorescence_ptr
  logical(kind=4), pointer :: ts_do_tf_iteration_ptr
  logical(kind=4), pointer :: ts_do_wladjusted_output_ptr
  logical(kind=4), pointer :: ts_do_toa_illumination_ptr
  logical(kind=4), pointer :: ts_do_boa_illumination_ptr
  logical(kind=4), dimension(:), pointer :: ts_do_albtrn_media_ptr
  logical(kind=4), pointer :: ts_do_planetary_problem_ptr
  logical(kind=4), pointer :: ts_do_mssts_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_do_fullrad_mode = .FALSE.
  ts_do_fullrad_mode_ptr => fortran_type_f%ts_do_fullrad_mode
  transfer_struct_c%ts_do_fullrad_mode = c_loc(ts_do_fullrad_mode_ptr)
  inquire(iolength=transfer_struct_c%ts_do_fullrad_mode_f_byte_size) fortran_type_f%ts_do_fullrad_mode
#ifdef ifort
  transfer_struct_c%ts_do_fullrad_mode_f_byte_size = transfer_struct_c%ts_do_fullrad_mode_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_thermal_emission = .FALSE.
  ts_do_thermal_emission_ptr => fortran_type_f%ts_do_thermal_emission
  transfer_struct_c%ts_do_thermal_emission = c_loc(ts_do_thermal_emission_ptr)
  inquire(iolength=transfer_struct_c%ts_do_thermal_emission_f_byte_size) fortran_type_f%ts_do_thermal_emission
#ifdef ifort
  transfer_struct_c%ts_do_thermal_emission_f_byte_size = transfer_struct_c%ts_do_thermal_emission_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_surface_emission = .FALSE.
  ts_do_surface_emission_ptr => fortran_type_f%ts_do_surface_emission
  transfer_struct_c%ts_do_surface_emission = c_loc(ts_do_surface_emission_ptr)
  inquire(iolength=transfer_struct_c%ts_do_surface_emission_f_byte_size) fortran_type_f%ts_do_surface_emission
#ifdef ifort
  transfer_struct_c%ts_do_surface_emission_f_byte_size = transfer_struct_c%ts_do_surface_emission_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_plane_parallel = .FALSE.
  ts_do_plane_parallel_ptr => fortran_type_f%ts_do_plane_parallel
  transfer_struct_c%ts_do_plane_parallel = c_loc(ts_do_plane_parallel_ptr)
  inquire(iolength=transfer_struct_c%ts_do_plane_parallel_f_byte_size) fortran_type_f%ts_do_plane_parallel
#ifdef ifort
  transfer_struct_c%ts_do_plane_parallel_f_byte_size = transfer_struct_c%ts_do_plane_parallel_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_brdf_surface = .FALSE.
  ts_do_brdf_surface_ptr => fortran_type_f%ts_do_brdf_surface
  transfer_struct_c%ts_do_brdf_surface = c_loc(ts_do_brdf_surface_ptr)
  inquire(iolength=transfer_struct_c%ts_do_brdf_surface_f_byte_size) fortran_type_f%ts_do_brdf_surface
#ifdef ifort
  transfer_struct_c%ts_do_brdf_surface_f_byte_size = transfer_struct_c%ts_do_brdf_surface_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_upwelling = .FALSE.
  ts_do_upwelling_ptr => fortran_type_f%ts_do_upwelling
  transfer_struct_c%ts_do_upwelling = c_loc(ts_do_upwelling_ptr)
  inquire(iolength=transfer_struct_c%ts_do_upwelling_f_byte_size) fortran_type_f%ts_do_upwelling
#ifdef ifort
  transfer_struct_c%ts_do_upwelling_f_byte_size = transfer_struct_c%ts_do_upwelling_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_dnwelling = .FALSE.
  ts_do_dnwelling_ptr => fortran_type_f%ts_do_dnwelling
  transfer_struct_c%ts_do_dnwelling = c_loc(ts_do_dnwelling_ptr)
  inquire(iolength=transfer_struct_c%ts_do_dnwelling_f_byte_size) fortran_type_f%ts_do_dnwelling
#ifdef ifort
  transfer_struct_c%ts_do_dnwelling_f_byte_size = transfer_struct_c%ts_do_dnwelling_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_toa_contribs = .FALSE.
  ts_do_toa_contribs_ptr => fortran_type_f%ts_do_toa_contribs
  transfer_struct_c%ts_do_toa_contribs = c_loc(ts_do_toa_contribs_ptr)
  inquire(iolength=transfer_struct_c%ts_do_toa_contribs_f_byte_size) fortran_type_f%ts_do_toa_contribs
#ifdef ifort
  transfer_struct_c%ts_do_toa_contribs_f_byte_size = transfer_struct_c%ts_do_toa_contribs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_surface_leaving = .FALSE.
  ts_do_surface_leaving_ptr => fortran_type_f%ts_do_surface_leaving
  transfer_struct_c%ts_do_surface_leaving = c_loc(ts_do_surface_leaving_ptr)
  inquire(iolength=transfer_struct_c%ts_do_surface_leaving_f_byte_size) fortran_type_f%ts_do_surface_leaving
#ifdef ifort
  transfer_struct_c%ts_do_surface_leaving_f_byte_size = transfer_struct_c%ts_do_surface_leaving_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sl_isotropic = .FALSE.
  ts_do_sl_isotropic_ptr => fortran_type_f%ts_do_sl_isotropic
  transfer_struct_c%ts_do_sl_isotropic = c_loc(ts_do_sl_isotropic_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sl_isotropic_f_byte_size) fortran_type_f%ts_do_sl_isotropic
#ifdef ifort
  transfer_struct_c%ts_do_sl_isotropic_f_byte_size = transfer_struct_c%ts_do_sl_isotropic_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_water_leaving = .FALSE.
  ts_do_water_leaving_ptr => fortran_type_f%ts_do_water_leaving
  transfer_struct_c%ts_do_water_leaving = c_loc(ts_do_water_leaving_ptr)
  inquire(iolength=transfer_struct_c%ts_do_water_leaving_f_byte_size) fortran_type_f%ts_do_water_leaving
#ifdef ifort
  transfer_struct_c%ts_do_water_leaving_f_byte_size = transfer_struct_c%ts_do_water_leaving_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_fluorescence = .FALSE.
  ts_do_fluorescence_ptr => fortran_type_f%ts_do_fluorescence
  transfer_struct_c%ts_do_fluorescence = c_loc(ts_do_fluorescence_ptr)
  inquire(iolength=transfer_struct_c%ts_do_fluorescence_f_byte_size) fortran_type_f%ts_do_fluorescence
#ifdef ifort
  transfer_struct_c%ts_do_fluorescence_f_byte_size = transfer_struct_c%ts_do_fluorescence_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_tf_iteration = .FALSE.
  ts_do_tf_iteration_ptr => fortran_type_f%ts_do_tf_iteration
  transfer_struct_c%ts_do_tf_iteration = c_loc(ts_do_tf_iteration_ptr)
  inquire(iolength=transfer_struct_c%ts_do_tf_iteration_f_byte_size) fortran_type_f%ts_do_tf_iteration
#ifdef ifort
  transfer_struct_c%ts_do_tf_iteration_f_byte_size = transfer_struct_c%ts_do_tf_iteration_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_wladjusted_output = .FALSE.
  ts_do_wladjusted_output_ptr => fortran_type_f%ts_do_wladjusted_output
  transfer_struct_c%ts_do_wladjusted_output = c_loc(ts_do_wladjusted_output_ptr)
  inquire(iolength=transfer_struct_c%ts_do_wladjusted_output_f_byte_size) fortran_type_f%ts_do_wladjusted_output
#ifdef ifort
  transfer_struct_c%ts_do_wladjusted_output_f_byte_size = transfer_struct_c%ts_do_wladjusted_output_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_toa_illumination = .FALSE.
  ts_do_toa_illumination_ptr => fortran_type_f%ts_do_toa_illumination
  transfer_struct_c%ts_do_toa_illumination = c_loc(ts_do_toa_illumination_ptr)
  inquire(iolength=transfer_struct_c%ts_do_toa_illumination_f_byte_size) fortran_type_f%ts_do_toa_illumination
#ifdef ifort
  transfer_struct_c%ts_do_toa_illumination_f_byte_size = transfer_struct_c%ts_do_toa_illumination_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_boa_illumination = .FALSE.
  ts_do_boa_illumination_ptr => fortran_type_f%ts_do_boa_illumination
  transfer_struct_c%ts_do_boa_illumination = c_loc(ts_do_boa_illumination_ptr)
  inquire(iolength=transfer_struct_c%ts_do_boa_illumination_f_byte_size) fortran_type_f%ts_do_boa_illumination
#ifdef ifort
  transfer_struct_c%ts_do_boa_illumination_f_byte_size = transfer_struct_c%ts_do_boa_illumination_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_albtrn_media = .FALSE.
  ts_do_albtrn_media_ptr => fortran_type_f%ts_do_albtrn_media
  transfer_struct_c%ts_do_albtrn_media = c_loc(ts_do_albtrn_media_ptr(&
    lbound(fortran_type_f%ts_do_albtrn_media,1)))
  inquire(iolength=transfer_struct_c%ts_do_albtrn_media_f_byte_size) fortran_type_f%ts_do_albtrn_media(&
    lbound(fortran_type_f%ts_do_albtrn_media,1))
#ifdef ifort
  transfer_struct_c%ts_do_albtrn_media_f_byte_size = transfer_struct_c%ts_do_albtrn_media_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_do_albtrn_media_f_shapes(1) = size(fortran_type_f%ts_do_albtrn_media, 1)
  
  
  fortran_type_f%ts_do_planetary_problem = .FALSE.
  ts_do_planetary_problem_ptr => fortran_type_f%ts_do_planetary_problem
  transfer_struct_c%ts_do_planetary_problem = c_loc(ts_do_planetary_problem_ptr)
  inquire(iolength=transfer_struct_c%ts_do_planetary_problem_f_byte_size) fortran_type_f%ts_do_planetary_problem
#ifdef ifort
  transfer_struct_c%ts_do_planetary_problem_f_byte_size = transfer_struct_c%ts_do_planetary_problem_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_mssts = .FALSE.
  ts_do_mssts_ptr => fortran_type_f%ts_do_mssts
  transfer_struct_c%ts_do_mssts = c_loc(ts_do_mssts_ptr)
  inquire(iolength=transfer_struct_c%ts_do_mssts_f_byte_size) fortran_type_f%ts_do_mssts
#ifdef ifort
  transfer_struct_c%ts_do_mssts_f_byte_size = transfer_struct_c%ts_do_mssts_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_boolean_c_init_only

subroutine lidort_fixed_boolean_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_boolean

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_boolean_c_destroy

subroutine lidort_fixed_boolean_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_boolean

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_boolean), pointer :: fortran_type_f_from
  type(lidort_fixed_boolean), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_fullrad_mode = fortran_type_f_from%ts_do_fullrad_mode
  fortran_type_f_to%ts_do_thermal_emission = fortran_type_f_from%ts_do_thermal_emission
  fortran_type_f_to%ts_do_surface_emission = fortran_type_f_from%ts_do_surface_emission
  fortran_type_f_to%ts_do_plane_parallel = fortran_type_f_from%ts_do_plane_parallel
  fortran_type_f_to%ts_do_brdf_surface = fortran_type_f_from%ts_do_brdf_surface
  fortran_type_f_to%ts_do_upwelling = fortran_type_f_from%ts_do_upwelling
  fortran_type_f_to%ts_do_dnwelling = fortran_type_f_from%ts_do_dnwelling
  fortran_type_f_to%ts_do_toa_contribs = fortran_type_f_from%ts_do_toa_contribs
  fortran_type_f_to%ts_do_surface_leaving = fortran_type_f_from%ts_do_surface_leaving
  fortran_type_f_to%ts_do_sl_isotropic = fortran_type_f_from%ts_do_sl_isotropic
  fortran_type_f_to%ts_do_water_leaving = fortran_type_f_from%ts_do_water_leaving
  fortran_type_f_to%ts_do_fluorescence = fortran_type_f_from%ts_do_fluorescence
  fortran_type_f_to%ts_do_tf_iteration = fortran_type_f_from%ts_do_tf_iteration
  fortran_type_f_to%ts_do_wladjusted_output = fortran_type_f_from%ts_do_wladjusted_output
  fortran_type_f_to%ts_do_toa_illumination = fortran_type_f_from%ts_do_toa_illumination
  fortran_type_f_to%ts_do_boa_illumination = fortran_type_f_from%ts_do_boa_illumination
  fortran_type_f_to%ts_do_albtrn_media = fortran_type_f_from%ts_do_albtrn_media
  fortran_type_f_to%ts_do_planetary_problem = fortran_type_f_from%ts_do_planetary_problem
  fortran_type_f_to%ts_do_mssts = fortran_type_f_from%ts_do_mssts
  

end subroutine lidort_fixed_boolean_c_copy

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_control_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_control

  type(lidort_fixed_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_control_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_control_c_alloc_init

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_control_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_taylor_order_ptr
  integer(c_int), pointer :: ts_nstreams_ptr
  integer(c_int), pointer :: ts_nlayers_ptr
  integer(c_int), pointer :: ts_nfinelayers_ptr
  integer(c_int), pointer :: ts_n_thermal_coeffs_ptr
  real(c_double), pointer :: ts_lidort_accuracy_ptr
  real(c_double), pointer :: ts_asymtx_tolerance_ptr
  integer(c_int), pointer :: ts_tf_maxiter_ptr
  real(c_double), pointer :: ts_tf_criterion_ptr
  real(c_double), pointer :: ts_toa_illumination_ptr
  real(c_double), pointer :: ts_boa_illumination_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_taylor_order = 0
  ts_taylor_order_ptr => fortran_type_f%ts_taylor_order
  transfer_struct_c%ts_taylor_order = c_loc(ts_taylor_order_ptr)
  inquire(iolength=transfer_struct_c%ts_taylor_order_f_byte_size) fortran_type_f%ts_taylor_order
#ifdef ifort
  transfer_struct_c%ts_taylor_order_f_byte_size = transfer_struct_c%ts_taylor_order_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_nstreams = 0
  ts_nstreams_ptr => fortran_type_f%ts_nstreams
  transfer_struct_c%ts_nstreams = c_loc(ts_nstreams_ptr)
  inquire(iolength=transfer_struct_c%ts_nstreams_f_byte_size) fortran_type_f%ts_nstreams
#ifdef ifort
  transfer_struct_c%ts_nstreams_f_byte_size = transfer_struct_c%ts_nstreams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_nlayers = 0
  ts_nlayers_ptr => fortran_type_f%ts_nlayers
  transfer_struct_c%ts_nlayers = c_loc(ts_nlayers_ptr)
  inquire(iolength=transfer_struct_c%ts_nlayers_f_byte_size) fortran_type_f%ts_nlayers
#ifdef ifort
  transfer_struct_c%ts_nlayers_f_byte_size = transfer_struct_c%ts_nlayers_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_nfinelayers = 0
  ts_nfinelayers_ptr => fortran_type_f%ts_nfinelayers
  transfer_struct_c%ts_nfinelayers = c_loc(ts_nfinelayers_ptr)
  inquire(iolength=transfer_struct_c%ts_nfinelayers_f_byte_size) fortran_type_f%ts_nfinelayers
#ifdef ifort
  transfer_struct_c%ts_nfinelayers_f_byte_size = transfer_struct_c%ts_nfinelayers_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_n_thermal_coeffs = 0
  ts_n_thermal_coeffs_ptr => fortran_type_f%ts_n_thermal_coeffs
  transfer_struct_c%ts_n_thermal_coeffs = c_loc(ts_n_thermal_coeffs_ptr)
  inquire(iolength=transfer_struct_c%ts_n_thermal_coeffs_f_byte_size) fortran_type_f%ts_n_thermal_coeffs
#ifdef ifort
  transfer_struct_c%ts_n_thermal_coeffs_f_byte_size = transfer_struct_c%ts_n_thermal_coeffs_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_lidort_accuracy = 0_fpk
  ts_lidort_accuracy_ptr => fortran_type_f%ts_lidort_accuracy
  transfer_struct_c%ts_lidort_accuracy = c_loc(ts_lidort_accuracy_ptr)
  inquire(iolength=transfer_struct_c%ts_lidort_accuracy_f_byte_size) fortran_type_f%ts_lidort_accuracy
#ifdef ifort
  transfer_struct_c%ts_lidort_accuracy_f_byte_size = transfer_struct_c%ts_lidort_accuracy_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_asymtx_tolerance = 0_fpk
  ts_asymtx_tolerance_ptr => fortran_type_f%ts_asymtx_tolerance
  transfer_struct_c%ts_asymtx_tolerance = c_loc(ts_asymtx_tolerance_ptr)
  inquire(iolength=transfer_struct_c%ts_asymtx_tolerance_f_byte_size) fortran_type_f%ts_asymtx_tolerance
#ifdef ifort
  transfer_struct_c%ts_asymtx_tolerance_f_byte_size = transfer_struct_c%ts_asymtx_tolerance_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_tf_maxiter = 0
  ts_tf_maxiter_ptr => fortran_type_f%ts_tf_maxiter
  transfer_struct_c%ts_tf_maxiter = c_loc(ts_tf_maxiter_ptr)
  inquire(iolength=transfer_struct_c%ts_tf_maxiter_f_byte_size) fortran_type_f%ts_tf_maxiter
#ifdef ifort
  transfer_struct_c%ts_tf_maxiter_f_byte_size = transfer_struct_c%ts_tf_maxiter_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_tf_criterion = 0_fpk
  ts_tf_criterion_ptr => fortran_type_f%ts_tf_criterion
  transfer_struct_c%ts_tf_criterion = c_loc(ts_tf_criterion_ptr)
  inquire(iolength=transfer_struct_c%ts_tf_criterion_f_byte_size) fortran_type_f%ts_tf_criterion
#ifdef ifort
  transfer_struct_c%ts_tf_criterion_f_byte_size = transfer_struct_c%ts_tf_criterion_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_toa_illumination = 0_fpk
  ts_toa_illumination_ptr => fortran_type_f%ts_toa_illumination
  transfer_struct_c%ts_toa_illumination = c_loc(ts_toa_illumination_ptr)
  inquire(iolength=transfer_struct_c%ts_toa_illumination_f_byte_size) fortran_type_f%ts_toa_illumination
#ifdef ifort
  transfer_struct_c%ts_toa_illumination_f_byte_size = transfer_struct_c%ts_toa_illumination_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_boa_illumination = 0_fpk
  ts_boa_illumination_ptr => fortran_type_f%ts_boa_illumination
  transfer_struct_c%ts_boa_illumination = c_loc(ts_boa_illumination_ptr)
  inquire(iolength=transfer_struct_c%ts_boa_illumination_f_byte_size) fortran_type_f%ts_boa_illumination
#ifdef ifort
  transfer_struct_c%ts_boa_illumination_f_byte_size = transfer_struct_c%ts_boa_illumination_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_control_c_init_only

subroutine lidort_fixed_control_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_control

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_control_c_destroy

subroutine lidort_fixed_control_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_control

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_control), pointer :: fortran_type_f_from
  type(lidort_fixed_control), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_taylor_order = fortran_type_f_from%ts_taylor_order
  fortran_type_f_to%ts_nstreams = fortran_type_f_from%ts_nstreams
  fortran_type_f_to%ts_nlayers = fortran_type_f_from%ts_nlayers
  fortran_type_f_to%ts_nfinelayers = fortran_type_f_from%ts_nfinelayers
  fortran_type_f_to%ts_n_thermal_coeffs = fortran_type_f_from%ts_n_thermal_coeffs
  fortran_type_f_to%ts_lidort_accuracy = fortran_type_f_from%ts_lidort_accuracy
  fortran_type_f_to%ts_asymtx_tolerance = fortran_type_f_from%ts_asymtx_tolerance
  fortran_type_f_to%ts_tf_maxiter = fortran_type_f_from%ts_tf_maxiter
  fortran_type_f_to%ts_tf_criterion = fortran_type_f_from%ts_tf_criterion
  fortran_type_f_to%ts_toa_illumination = fortran_type_f_from%ts_toa_illumination
  fortran_type_f_to%ts_boa_illumination = fortran_type_f_from%ts_boa_illumination
  

end subroutine lidort_fixed_control_c_copy

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_sunrays_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_sunrays

  type(lidort_fixed_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_sunrays_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_sunrays_c_alloc_init

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_sunrays_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  real(c_double), pointer :: ts_flux_factor_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_flux_factor = 0_fpk
  ts_flux_factor_ptr => fortran_type_f%ts_flux_factor
  transfer_struct_c%ts_flux_factor = c_loc(ts_flux_factor_ptr)
  inquire(iolength=transfer_struct_c%ts_flux_factor_f_byte_size) fortran_type_f%ts_flux_factor
#ifdef ifort
  transfer_struct_c%ts_flux_factor_f_byte_size = transfer_struct_c%ts_flux_factor_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_sunrays_c_init_only

subroutine lidort_fixed_sunrays_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_sunrays_c_destroy

subroutine lidort_fixed_sunrays_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_sunrays), pointer :: fortran_type_f_from
  type(lidort_fixed_sunrays), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_flux_factor = fortran_type_f_from%ts_flux_factor
  

end subroutine lidort_fixed_sunrays_c_copy

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_uservalues_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_uservalues

  type(lidort_fixed_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_uservalues_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_uservalues_c_alloc_init

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_uservalues_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_n_user_levels_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_n_user_levels = 0
  ts_n_user_levels_ptr => fortran_type_f%ts_n_user_levels
  transfer_struct_c%ts_n_user_levels = c_loc(ts_n_user_levels_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_levels_f_byte_size) fortran_type_f%ts_n_user_levels
#ifdef ifort
  transfer_struct_c%ts_n_user_levels_f_byte_size = transfer_struct_c%ts_n_user_levels_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_uservalues_c_init_only

subroutine lidort_fixed_uservalues_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_uservalues_c_destroy

subroutine lidort_fixed_uservalues_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_uservalues), pointer :: fortran_type_f_from
  type(lidort_fixed_uservalues), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_n_user_levels = fortran_type_f_from%ts_n_user_levels
  

end subroutine lidort_fixed_uservalues_c_copy

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_chapman_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_chapman

  type(lidort_fixed_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_chapman_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_chapman_c_alloc_init

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_chapman_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  real(c_double), dimension(:), pointer :: ts_height_grid_ptr
  real(c_double), dimension(:), pointer :: ts_pressure_grid_ptr
  real(c_double), dimension(:), pointer :: ts_temperature_grid_ptr
  integer(c_int), dimension(:), pointer :: ts_finegrid_ptr
  real(c_double), pointer :: ts_rfindex_parameter_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_height_grid = 0_fpk
  ts_height_grid_ptr => fortran_type_f%ts_height_grid
  transfer_struct_c%ts_height_grid = c_loc(ts_height_grid_ptr(&
    lbound(fortran_type_f%ts_height_grid,1)))
  inquire(iolength=transfer_struct_c%ts_height_grid_f_byte_size) fortran_type_f%ts_height_grid(&
    lbound(fortran_type_f%ts_height_grid,1))
#ifdef ifort
  transfer_struct_c%ts_height_grid_f_byte_size = transfer_struct_c%ts_height_grid_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_height_grid_f_shapes(1) = size(fortran_type_f%ts_height_grid, 1)
  
  
  fortran_type_f%ts_pressure_grid = 0_fpk
  ts_pressure_grid_ptr => fortran_type_f%ts_pressure_grid
  transfer_struct_c%ts_pressure_grid = c_loc(ts_pressure_grid_ptr(&
    lbound(fortran_type_f%ts_pressure_grid,1)))
  inquire(iolength=transfer_struct_c%ts_pressure_grid_f_byte_size) fortran_type_f%ts_pressure_grid(&
    lbound(fortran_type_f%ts_pressure_grid,1))
#ifdef ifort
  transfer_struct_c%ts_pressure_grid_f_byte_size = transfer_struct_c%ts_pressure_grid_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_pressure_grid_f_shapes(1) = size(fortran_type_f%ts_pressure_grid, 1)
  
  
  fortran_type_f%ts_temperature_grid = 0_fpk
  ts_temperature_grid_ptr => fortran_type_f%ts_temperature_grid
  transfer_struct_c%ts_temperature_grid = c_loc(ts_temperature_grid_ptr(&
    lbound(fortran_type_f%ts_temperature_grid,1)))
  inquire(iolength=transfer_struct_c%ts_temperature_grid_f_byte_size) fortran_type_f%ts_temperature_grid(&
    lbound(fortran_type_f%ts_temperature_grid,1))
#ifdef ifort
  transfer_struct_c%ts_temperature_grid_f_byte_size = transfer_struct_c%ts_temperature_grid_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_temperature_grid_f_shapes(1) = size(fortran_type_f%ts_temperature_grid, 1)
  
  
  fortran_type_f%ts_finegrid = 0
  ts_finegrid_ptr => fortran_type_f%ts_finegrid
  transfer_struct_c%ts_finegrid = c_loc(ts_finegrid_ptr(&
    lbound(fortran_type_f%ts_finegrid,1)))
  inquire(iolength=transfer_struct_c%ts_finegrid_f_byte_size) fortran_type_f%ts_finegrid(&
    lbound(fortran_type_f%ts_finegrid,1))
#ifdef ifort
  transfer_struct_c%ts_finegrid_f_byte_size = transfer_struct_c%ts_finegrid_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_finegrid_f_shapes(1) = size(fortran_type_f%ts_finegrid, 1)
  
  
  fortran_type_f%ts_rfindex_parameter = 0_fpk
  ts_rfindex_parameter_ptr => fortran_type_f%ts_rfindex_parameter
  transfer_struct_c%ts_rfindex_parameter = c_loc(ts_rfindex_parameter_ptr)
  inquire(iolength=transfer_struct_c%ts_rfindex_parameter_f_byte_size) fortran_type_f%ts_rfindex_parameter
#ifdef ifort
  transfer_struct_c%ts_rfindex_parameter_f_byte_size = transfer_struct_c%ts_rfindex_parameter_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_chapman_c_init_only

subroutine lidort_fixed_chapman_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_chapman

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_chapman_c_destroy

subroutine lidort_fixed_chapman_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_chapman

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_chapman), pointer :: fortran_type_f_from
  type(lidort_fixed_chapman), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_height_grid = fortran_type_f_from%ts_height_grid
  fortran_type_f_to%ts_pressure_grid = fortran_type_f_from%ts_pressure_grid
  fortran_type_f_to%ts_temperature_grid = fortran_type_f_from%ts_temperature_grid
  fortran_type_f_to%ts_finegrid = fortran_type_f_from%ts_finegrid
  fortran_type_f_to%ts_rfindex_parameter = fortran_type_f_from%ts_rfindex_parameter
  

end subroutine lidort_fixed_chapman_c_copy

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_optical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_optical

  type(lidort_fixed_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_optical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_optical_c_alloc_init

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_optical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  real(c_double), dimension(:), pointer :: ts_deltau_vert_input_ptr
  real(c_double), dimension(:,:), pointer :: ts_phasmoms_total_input_ptr
  real(c_double), dimension(:,:), pointer :: ts_phasfunc_input_up_ptr
  real(c_double), dimension(:,:), pointer :: ts_phasfunc_input_dn_ptr
  real(c_double), pointer :: ts_lambertian_albedo_ptr
  real(c_double), dimension(:), pointer :: ts_thermal_bb_input_ptr
  real(c_double), pointer :: ts_surface_bb_input_ptr
  real(c_double), pointer :: ts_atmos_wavelength_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_deltau_vert_input = 0_fpk
  ts_deltau_vert_input_ptr => fortran_type_f%ts_deltau_vert_input
  transfer_struct_c%ts_deltau_vert_input = c_loc(ts_deltau_vert_input_ptr(&
    lbound(fortran_type_f%ts_deltau_vert_input,1)))
  inquire(iolength=transfer_struct_c%ts_deltau_vert_input_f_byte_size) fortran_type_f%ts_deltau_vert_input(&
    lbound(fortran_type_f%ts_deltau_vert_input,1))
#ifdef ifort
  transfer_struct_c%ts_deltau_vert_input_f_byte_size = transfer_struct_c%ts_deltau_vert_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_deltau_vert_input_f_shapes(1) = size(fortran_type_f%ts_deltau_vert_input, 1)
  
  
  fortran_type_f%ts_phasmoms_total_input = 0_fpk
  ts_phasmoms_total_input_ptr => fortran_type_f%ts_phasmoms_total_input
  transfer_struct_c%ts_phasmoms_total_input = c_loc(ts_phasmoms_total_input_ptr(&
    lbound(fortran_type_f%ts_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_phasmoms_total_input,2)))
  inquire(iolength=transfer_struct_c%ts_phasmoms_total_input_f_byte_size) fortran_type_f%ts_phasmoms_total_input(&
    lbound(fortran_type_f%ts_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_phasmoms_total_input,2))
#ifdef ifort
  transfer_struct_c%ts_phasmoms_total_input_f_byte_size = transfer_struct_c%ts_phasmoms_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_phasmoms_total_input_f_shapes(1) = size(fortran_type_f%ts_phasmoms_total_input, 1)
  transfer_struct_c%ts_phasmoms_total_input_f_shapes(2) = size(fortran_type_f%ts_phasmoms_total_input, 2)
  
  
  fortran_type_f%ts_phasfunc_input_up = 0_fpk
  ts_phasfunc_input_up_ptr => fortran_type_f%ts_phasfunc_input_up
  transfer_struct_c%ts_phasfunc_input_up = c_loc(ts_phasfunc_input_up_ptr(&
    lbound(fortran_type_f%ts_phasfunc_input_up,1),&
    lbound(fortran_type_f%ts_phasfunc_input_up,2)))
  inquire(iolength=transfer_struct_c%ts_phasfunc_input_up_f_byte_size) fortran_type_f%ts_phasfunc_input_up(&
    lbound(fortran_type_f%ts_phasfunc_input_up,1),&
    lbound(fortran_type_f%ts_phasfunc_input_up,2))
#ifdef ifort
  transfer_struct_c%ts_phasfunc_input_up_f_byte_size = transfer_struct_c%ts_phasfunc_input_up_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_phasfunc_input_up_f_shapes(1) = size(fortran_type_f%ts_phasfunc_input_up, 1)
  transfer_struct_c%ts_phasfunc_input_up_f_shapes(2) = size(fortran_type_f%ts_phasfunc_input_up, 2)
  
  
  fortran_type_f%ts_phasfunc_input_dn = 0_fpk
  ts_phasfunc_input_dn_ptr => fortran_type_f%ts_phasfunc_input_dn
  transfer_struct_c%ts_phasfunc_input_dn = c_loc(ts_phasfunc_input_dn_ptr(&
    lbound(fortran_type_f%ts_phasfunc_input_dn,1),&
    lbound(fortran_type_f%ts_phasfunc_input_dn,2)))
  inquire(iolength=transfer_struct_c%ts_phasfunc_input_dn_f_byte_size) fortran_type_f%ts_phasfunc_input_dn(&
    lbound(fortran_type_f%ts_phasfunc_input_dn,1),&
    lbound(fortran_type_f%ts_phasfunc_input_dn,2))
#ifdef ifort
  transfer_struct_c%ts_phasfunc_input_dn_f_byte_size = transfer_struct_c%ts_phasfunc_input_dn_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_phasfunc_input_dn_f_shapes(1) = size(fortran_type_f%ts_phasfunc_input_dn, 1)
  transfer_struct_c%ts_phasfunc_input_dn_f_shapes(2) = size(fortran_type_f%ts_phasfunc_input_dn, 2)
  
  
  fortran_type_f%ts_lambertian_albedo = 0_fpk
  ts_lambertian_albedo_ptr => fortran_type_f%ts_lambertian_albedo
  transfer_struct_c%ts_lambertian_albedo = c_loc(ts_lambertian_albedo_ptr)
  inquire(iolength=transfer_struct_c%ts_lambertian_albedo_f_byte_size) fortran_type_f%ts_lambertian_albedo
#ifdef ifort
  transfer_struct_c%ts_lambertian_albedo_f_byte_size = transfer_struct_c%ts_lambertian_albedo_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_thermal_bb_input = 0_fpk
  ts_thermal_bb_input_ptr => fortran_type_f%ts_thermal_bb_input
  transfer_struct_c%ts_thermal_bb_input = c_loc(ts_thermal_bb_input_ptr(&
    lbound(fortran_type_f%ts_thermal_bb_input,1)))
  inquire(iolength=transfer_struct_c%ts_thermal_bb_input_f_byte_size) fortran_type_f%ts_thermal_bb_input(&
    lbound(fortran_type_f%ts_thermal_bb_input,1))
#ifdef ifort
  transfer_struct_c%ts_thermal_bb_input_f_byte_size = transfer_struct_c%ts_thermal_bb_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_thermal_bb_input_f_shapes(1) = size(fortran_type_f%ts_thermal_bb_input, 1)
  
  
  fortran_type_f%ts_surface_bb_input = 0_fpk
  ts_surface_bb_input_ptr => fortran_type_f%ts_surface_bb_input
  transfer_struct_c%ts_surface_bb_input = c_loc(ts_surface_bb_input_ptr)
  inquire(iolength=transfer_struct_c%ts_surface_bb_input_f_byte_size) fortran_type_f%ts_surface_bb_input
#ifdef ifort
  transfer_struct_c%ts_surface_bb_input_f_byte_size = transfer_struct_c%ts_surface_bb_input_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_atmos_wavelength = 0_fpk
  ts_atmos_wavelength_ptr => fortran_type_f%ts_atmos_wavelength
  transfer_struct_c%ts_atmos_wavelength = c_loc(ts_atmos_wavelength_ptr)
  inquire(iolength=transfer_struct_c%ts_atmos_wavelength_f_byte_size) fortran_type_f%ts_atmos_wavelength
#ifdef ifort
  transfer_struct_c%ts_atmos_wavelength_f_byte_size = transfer_struct_c%ts_atmos_wavelength_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_optical_c_init_only

subroutine lidort_fixed_optical_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_optical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_optical_c_destroy

subroutine lidort_fixed_optical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_optical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_optical), pointer :: fortran_type_f_from
  type(lidort_fixed_optical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_deltau_vert_input = fortran_type_f_from%ts_deltau_vert_input
  fortran_type_f_to%ts_phasmoms_total_input = fortran_type_f_from%ts_phasmoms_total_input
  fortran_type_f_to%ts_phasfunc_input_up = fortran_type_f_from%ts_phasfunc_input_up
  fortran_type_f_to%ts_phasfunc_input_dn = fortran_type_f_from%ts_phasfunc_input_dn
  fortran_type_f_to%ts_lambertian_albedo = fortran_type_f_from%ts_lambertian_albedo
  fortran_type_f_to%ts_thermal_bb_input = fortran_type_f_from%ts_thermal_bb_input
  fortran_type_f_to%ts_surface_bb_input = fortran_type_f_from%ts_surface_bb_input
  fortran_type_f_to%ts_atmos_wavelength = fortran_type_f_from%ts_atmos_wavelength
  

end subroutine lidort_fixed_optical_c_copy

! Links to type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_write_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(lidort_fixed_write_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_write), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_write_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_write_c_alloc_init

! Links to type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_write_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_write_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_write), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_debug_write_ptr
  logical(kind=4), pointer :: ts_do_write_input_ptr
  logical(kind=4), pointer :: ts_do_write_scenario_ptr
  logical(kind=4), pointer :: ts_do_write_fourier_ptr
  logical(kind=4), pointer :: ts_do_write_results_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_do_debug_write = .FALSE.
  ts_do_debug_write_ptr => fortran_type_f%ts_do_debug_write
  transfer_struct_c%ts_do_debug_write = c_loc(ts_do_debug_write_ptr)
  inquire(iolength=transfer_struct_c%ts_do_debug_write_f_byte_size) fortran_type_f%ts_do_debug_write
#ifdef ifort
  transfer_struct_c%ts_do_debug_write_f_byte_size = transfer_struct_c%ts_do_debug_write_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_write_input = .FALSE.
  ts_do_write_input_ptr => fortran_type_f%ts_do_write_input
  transfer_struct_c%ts_do_write_input = c_loc(ts_do_write_input_ptr)
  inquire(iolength=transfer_struct_c%ts_do_write_input_f_byte_size) fortran_type_f%ts_do_write_input
#ifdef ifort
  transfer_struct_c%ts_do_write_input_f_byte_size = transfer_struct_c%ts_do_write_input_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_input_write_filename = ''
  transfer_struct_c%ts_input_write_filename_f_len = len(fortran_type_f%ts_input_write_filename)
  
  
  fortran_type_f%ts_do_write_scenario = .FALSE.
  ts_do_write_scenario_ptr => fortran_type_f%ts_do_write_scenario
  transfer_struct_c%ts_do_write_scenario = c_loc(ts_do_write_scenario_ptr)
  inquire(iolength=transfer_struct_c%ts_do_write_scenario_f_byte_size) fortran_type_f%ts_do_write_scenario
#ifdef ifort
  transfer_struct_c%ts_do_write_scenario_f_byte_size = transfer_struct_c%ts_do_write_scenario_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_scenario_write_filename = ''
  transfer_struct_c%ts_scenario_write_filename_f_len = len(fortran_type_f%ts_scenario_write_filename)
  
  
  fortran_type_f%ts_do_write_fourier = .FALSE.
  ts_do_write_fourier_ptr => fortran_type_f%ts_do_write_fourier
  transfer_struct_c%ts_do_write_fourier = c_loc(ts_do_write_fourier_ptr)
  inquire(iolength=transfer_struct_c%ts_do_write_fourier_f_byte_size) fortran_type_f%ts_do_write_fourier
#ifdef ifort
  transfer_struct_c%ts_do_write_fourier_f_byte_size = transfer_struct_c%ts_do_write_fourier_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_fourier_write_filename = ''
  transfer_struct_c%ts_fourier_write_filename_f_len = len(fortran_type_f%ts_fourier_write_filename)
  
  
  fortran_type_f%ts_do_write_results = .FALSE.
  ts_do_write_results_ptr => fortran_type_f%ts_do_write_results
  transfer_struct_c%ts_do_write_results = c_loc(ts_do_write_results_ptr)
  inquire(iolength=transfer_struct_c%ts_do_write_results_f_byte_size) fortran_type_f%ts_do_write_results
#ifdef ifort
  transfer_struct_c%ts_do_write_results_f_byte_size = transfer_struct_c%ts_do_write_results_f_byte_size * 4
#endif
  
  
  fortran_type_f%ts_results_write_filename = ''
  transfer_struct_c%ts_results_write_filename_f_len = len(fortran_type_f%ts_results_write_filename)
  
  
end subroutine lidort_fixed_write_c_init_only

subroutine lidort_fixed_write_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_write), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_write_c_destroy

subroutine lidort_fixed_write_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_write), pointer :: fortran_type_f_from
  type(lidort_fixed_write), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_debug_write = fortran_type_f_from%ts_do_debug_write
  fortran_type_f_to%ts_do_write_input = fortran_type_f_from%ts_do_write_input
  fortran_type_f_to%ts_input_write_filename = fortran_type_f_from%ts_input_write_filename
  fortran_type_f_to%ts_do_write_scenario = fortran_type_f_from%ts_do_write_scenario
  fortran_type_f_to%ts_scenario_write_filename = fortran_type_f_from%ts_scenario_write_filename
  fortran_type_f_to%ts_do_write_fourier = fortran_type_f_from%ts_do_write_fourier
  fortran_type_f_to%ts_fourier_write_filename = fortran_type_f_from%ts_fourier_write_filename
  fortran_type_f_to%ts_do_write_results = fortran_type_f_from%ts_do_write_results
  fortran_type_f_to%ts_results_write_filename = fortran_type_f_from%ts_results_write_filename
  

end subroutine lidort_fixed_write_c_copy

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_inputs

  type(lidort_fixed_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_inputs_c_alloc_init

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_fixed_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  type(lidort_fixed_boolean), pointer :: bool_ptr
  type(lidort_fixed_control), pointer :: cont_ptr
  type(lidort_fixed_sunrays), pointer :: sunrays_ptr
  type(lidort_fixed_uservalues), pointer :: userval_ptr
  type(lidort_fixed_chapman), pointer :: chapman_ptr
  type(lidort_fixed_optical), pointer :: optical_ptr
  type(lidort_fixed_write), pointer :: write_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  bool_ptr => fortran_type_f%bool
  transfer_struct_c%f_bool = c_loc(bool_ptr)
  inquire(iolength=transfer_struct_c%bool_f_byte_size) fortran_type_f%bool
#ifdef ifort
  transfer_struct_c%bool_f_byte_size = transfer_struct_c%bool_f_byte_size * 4
#endif
  
  
  cont_ptr => fortran_type_f%cont
  transfer_struct_c%cont = c_loc(cont_ptr)
  inquire(iolength=transfer_struct_c%cont_f_byte_size) fortran_type_f%cont
#ifdef ifort
  transfer_struct_c%cont_f_byte_size = transfer_struct_c%cont_f_byte_size * 4
#endif
  
  
  sunrays_ptr => fortran_type_f%sunrays
  transfer_struct_c%sunrays = c_loc(sunrays_ptr)
  inquire(iolength=transfer_struct_c%sunrays_f_byte_size) fortran_type_f%sunrays
#ifdef ifort
  transfer_struct_c%sunrays_f_byte_size = transfer_struct_c%sunrays_f_byte_size * 4
#endif
  
  
  userval_ptr => fortran_type_f%userval
  transfer_struct_c%userval = c_loc(userval_ptr)
  inquire(iolength=transfer_struct_c%userval_f_byte_size) fortran_type_f%userval
#ifdef ifort
  transfer_struct_c%userval_f_byte_size = transfer_struct_c%userval_f_byte_size * 4
#endif
  
  
  chapman_ptr => fortran_type_f%chapman
  transfer_struct_c%chapman = c_loc(chapman_ptr)
  inquire(iolength=transfer_struct_c%chapman_f_byte_size) fortran_type_f%chapman
#ifdef ifort
  transfer_struct_c%chapman_f_byte_size = transfer_struct_c%chapman_f_byte_size * 4
#endif
  
  
  optical_ptr => fortran_type_f%optical
  transfer_struct_c%optical = c_loc(optical_ptr)
  inquire(iolength=transfer_struct_c%optical_f_byte_size) fortran_type_f%optical
#ifdef ifort
  transfer_struct_c%optical_f_byte_size = transfer_struct_c%optical_f_byte_size * 4
#endif
  
  
  write_ptr => fortran_type_f%write
  transfer_struct_c%write = c_loc(write_ptr)
  inquire(iolength=transfer_struct_c%write_f_byte_size) fortran_type_f%write
#ifdef ifort
  transfer_struct_c%write_f_byte_size = transfer_struct_c%write_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_inputs_c_init_only

subroutine lidort_fixed_inputs_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_inputs_c_destroy

subroutine lidort_fixed_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_inputs), pointer :: fortran_type_f_from
  type(lidort_fixed_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bool = fortran_type_f_from%bool
  fortran_type_f_to%cont = fortran_type_f_from%cont
  fortran_type_f_to%sunrays = fortran_type_f_from%sunrays
  fortran_type_f_to%userval = fortran_type_f_from%userval
  fortran_type_f_to%chapman = fortran_type_f_from%chapman
  fortran_type_f_to%optical = fortran_type_f_from%optical
  fortran_type_f_to%write = fortran_type_f_from%write
  

end subroutine lidort_fixed_inputs_c_copy

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_boolean_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_boolean

  type(lidort_modified_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_boolean_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_boolean_c_alloc_init

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_boolean_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_focorr_ptr
  logical(kind=4), pointer :: ts_do_focorr_external_ptr
  logical(kind=4), pointer :: ts_do_focorr_nadir_ptr
  logical(kind=4), pointer :: ts_do_focorr_outgoing_ptr
  logical(kind=4), pointer :: ts_do_sscorr_truncation_ptr
  logical(kind=4), pointer :: ts_do_sscorr_usephasfunc_ptr
  logical(kind=4), pointer :: ts_do_external_wleave_ptr
  logical(kind=4), pointer :: ts_do_double_convtest_ptr
  logical(kind=4), pointer :: ts_do_solar_sources_ptr
  logical(kind=4), pointer :: ts_do_refractive_geometry_ptr
  logical(kind=4), pointer :: ts_do_chapman_function_ptr
  logical(kind=4), pointer :: ts_do_rayleigh_only_ptr
  logical(kind=4), pointer :: ts_do_isotropic_only_ptr
  logical(kind=4), pointer :: ts_do_no_azimuth_ptr
  logical(kind=4), pointer :: ts_do_all_fourier_ptr
  logical(kind=4), pointer :: ts_do_deltam_scaling_ptr
  logical(kind=4), pointer :: ts_do_solution_saving_ptr
  logical(kind=4), pointer :: ts_do_bvp_telescoping_ptr
  logical(kind=4), pointer :: ts_do_user_streams_ptr
  logical(kind=4), pointer :: ts_do_additional_mvout_ptr
  logical(kind=4), pointer :: ts_do_mvout_only_ptr
  logical(kind=4), pointer :: ts_do_thermal_transonly_ptr
  logical(kind=4), pointer :: ts_do_observation_geometry_ptr
  logical(kind=4), pointer :: ts_do_doublet_geometry_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_do_focorr = .FALSE.
  ts_do_focorr_ptr => fortran_type_f%ts_do_focorr
  transfer_struct_c%ts_do_focorr = c_loc(ts_do_focorr_ptr)
  inquire(iolength=transfer_struct_c%ts_do_focorr_f_byte_size) fortran_type_f%ts_do_focorr
#ifdef ifort
  transfer_struct_c%ts_do_focorr_f_byte_size = transfer_struct_c%ts_do_focorr_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_focorr_external = .FALSE.
  ts_do_focorr_external_ptr => fortran_type_f%ts_do_focorr_external
  transfer_struct_c%ts_do_focorr_external = c_loc(ts_do_focorr_external_ptr)
  inquire(iolength=transfer_struct_c%ts_do_focorr_external_f_byte_size) fortran_type_f%ts_do_focorr_external
#ifdef ifort
  transfer_struct_c%ts_do_focorr_external_f_byte_size = transfer_struct_c%ts_do_focorr_external_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_focorr_nadir = .FALSE.
  ts_do_focorr_nadir_ptr => fortran_type_f%ts_do_focorr_nadir
  transfer_struct_c%ts_do_focorr_nadir = c_loc(ts_do_focorr_nadir_ptr)
  inquire(iolength=transfer_struct_c%ts_do_focorr_nadir_f_byte_size) fortran_type_f%ts_do_focorr_nadir
#ifdef ifort
  transfer_struct_c%ts_do_focorr_nadir_f_byte_size = transfer_struct_c%ts_do_focorr_nadir_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_focorr_outgoing = .FALSE.
  ts_do_focorr_outgoing_ptr => fortran_type_f%ts_do_focorr_outgoing
  transfer_struct_c%ts_do_focorr_outgoing = c_loc(ts_do_focorr_outgoing_ptr)
  inquire(iolength=transfer_struct_c%ts_do_focorr_outgoing_f_byte_size) fortran_type_f%ts_do_focorr_outgoing
#ifdef ifort
  transfer_struct_c%ts_do_focorr_outgoing_f_byte_size = transfer_struct_c%ts_do_focorr_outgoing_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sscorr_truncation = .FALSE.
  ts_do_sscorr_truncation_ptr => fortran_type_f%ts_do_sscorr_truncation
  transfer_struct_c%ts_do_sscorr_truncation = c_loc(ts_do_sscorr_truncation_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sscorr_truncation_f_byte_size) fortran_type_f%ts_do_sscorr_truncation
#ifdef ifort
  transfer_struct_c%ts_do_sscorr_truncation_f_byte_size = transfer_struct_c%ts_do_sscorr_truncation_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sscorr_usephasfunc = .FALSE.
  ts_do_sscorr_usephasfunc_ptr => fortran_type_f%ts_do_sscorr_usephasfunc
  transfer_struct_c%ts_do_sscorr_usephasfunc = c_loc(ts_do_sscorr_usephasfunc_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sscorr_usephasfunc_f_byte_size) fortran_type_f%ts_do_sscorr_usephasfunc
#ifdef ifort
  transfer_struct_c%ts_do_sscorr_usephasfunc_f_byte_size = transfer_struct_c%ts_do_sscorr_usephasfunc_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_external_wleave = .FALSE.
  ts_do_external_wleave_ptr => fortran_type_f%ts_do_external_wleave
  transfer_struct_c%ts_do_external_wleave = c_loc(ts_do_external_wleave_ptr)
  inquire(iolength=transfer_struct_c%ts_do_external_wleave_f_byte_size) fortran_type_f%ts_do_external_wleave
#ifdef ifort
  transfer_struct_c%ts_do_external_wleave_f_byte_size = transfer_struct_c%ts_do_external_wleave_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_double_convtest = .FALSE.
  ts_do_double_convtest_ptr => fortran_type_f%ts_do_double_convtest
  transfer_struct_c%ts_do_double_convtest = c_loc(ts_do_double_convtest_ptr)
  inquire(iolength=transfer_struct_c%ts_do_double_convtest_f_byte_size) fortran_type_f%ts_do_double_convtest
#ifdef ifort
  transfer_struct_c%ts_do_double_convtest_f_byte_size = transfer_struct_c%ts_do_double_convtest_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_solar_sources = .FALSE.
  ts_do_solar_sources_ptr => fortran_type_f%ts_do_solar_sources
  transfer_struct_c%ts_do_solar_sources = c_loc(ts_do_solar_sources_ptr)
  inquire(iolength=transfer_struct_c%ts_do_solar_sources_f_byte_size) fortran_type_f%ts_do_solar_sources
#ifdef ifort
  transfer_struct_c%ts_do_solar_sources_f_byte_size = transfer_struct_c%ts_do_solar_sources_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_refractive_geometry = .FALSE.
  ts_do_refractive_geometry_ptr => fortran_type_f%ts_do_refractive_geometry
  transfer_struct_c%ts_do_refractive_geometry = c_loc(ts_do_refractive_geometry_ptr)
  inquire(iolength=transfer_struct_c%ts_do_refractive_geometry_f_byte_size) fortran_type_f%ts_do_refractive_geometry
#ifdef ifort
  transfer_struct_c%ts_do_refractive_geometry_f_byte_size = transfer_struct_c%ts_do_refractive_geometry_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_chapman_function = .FALSE.
  ts_do_chapman_function_ptr => fortran_type_f%ts_do_chapman_function
  transfer_struct_c%ts_do_chapman_function = c_loc(ts_do_chapman_function_ptr)
  inquire(iolength=transfer_struct_c%ts_do_chapman_function_f_byte_size) fortran_type_f%ts_do_chapman_function
#ifdef ifort
  transfer_struct_c%ts_do_chapman_function_f_byte_size = transfer_struct_c%ts_do_chapman_function_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_rayleigh_only = .FALSE.
  ts_do_rayleigh_only_ptr => fortran_type_f%ts_do_rayleigh_only
  transfer_struct_c%ts_do_rayleigh_only = c_loc(ts_do_rayleigh_only_ptr)
  inquire(iolength=transfer_struct_c%ts_do_rayleigh_only_f_byte_size) fortran_type_f%ts_do_rayleigh_only
#ifdef ifort
  transfer_struct_c%ts_do_rayleigh_only_f_byte_size = transfer_struct_c%ts_do_rayleigh_only_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_isotropic_only = .FALSE.
  ts_do_isotropic_only_ptr => fortran_type_f%ts_do_isotropic_only
  transfer_struct_c%ts_do_isotropic_only = c_loc(ts_do_isotropic_only_ptr)
  inquire(iolength=transfer_struct_c%ts_do_isotropic_only_f_byte_size) fortran_type_f%ts_do_isotropic_only
#ifdef ifort
  transfer_struct_c%ts_do_isotropic_only_f_byte_size = transfer_struct_c%ts_do_isotropic_only_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_no_azimuth = .FALSE.
  ts_do_no_azimuth_ptr => fortran_type_f%ts_do_no_azimuth
  transfer_struct_c%ts_do_no_azimuth = c_loc(ts_do_no_azimuth_ptr)
  inquire(iolength=transfer_struct_c%ts_do_no_azimuth_f_byte_size) fortran_type_f%ts_do_no_azimuth
#ifdef ifort
  transfer_struct_c%ts_do_no_azimuth_f_byte_size = transfer_struct_c%ts_do_no_azimuth_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_all_fourier = .FALSE.
  ts_do_all_fourier_ptr => fortran_type_f%ts_do_all_fourier
  transfer_struct_c%ts_do_all_fourier = c_loc(ts_do_all_fourier_ptr)
  inquire(iolength=transfer_struct_c%ts_do_all_fourier_f_byte_size) fortran_type_f%ts_do_all_fourier
#ifdef ifort
  transfer_struct_c%ts_do_all_fourier_f_byte_size = transfer_struct_c%ts_do_all_fourier_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_deltam_scaling = .FALSE.
  ts_do_deltam_scaling_ptr => fortran_type_f%ts_do_deltam_scaling
  transfer_struct_c%ts_do_deltam_scaling = c_loc(ts_do_deltam_scaling_ptr)
  inquire(iolength=transfer_struct_c%ts_do_deltam_scaling_f_byte_size) fortran_type_f%ts_do_deltam_scaling
#ifdef ifort
  transfer_struct_c%ts_do_deltam_scaling_f_byte_size = transfer_struct_c%ts_do_deltam_scaling_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_solution_saving = .FALSE.
  ts_do_solution_saving_ptr => fortran_type_f%ts_do_solution_saving
  transfer_struct_c%ts_do_solution_saving = c_loc(ts_do_solution_saving_ptr)
  inquire(iolength=transfer_struct_c%ts_do_solution_saving_f_byte_size) fortran_type_f%ts_do_solution_saving
#ifdef ifort
  transfer_struct_c%ts_do_solution_saving_f_byte_size = transfer_struct_c%ts_do_solution_saving_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_bvp_telescoping = .FALSE.
  ts_do_bvp_telescoping_ptr => fortran_type_f%ts_do_bvp_telescoping
  transfer_struct_c%ts_do_bvp_telescoping = c_loc(ts_do_bvp_telescoping_ptr)
  inquire(iolength=transfer_struct_c%ts_do_bvp_telescoping_f_byte_size) fortran_type_f%ts_do_bvp_telescoping
#ifdef ifort
  transfer_struct_c%ts_do_bvp_telescoping_f_byte_size = transfer_struct_c%ts_do_bvp_telescoping_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_user_streams = .FALSE.
  ts_do_user_streams_ptr => fortran_type_f%ts_do_user_streams
  transfer_struct_c%ts_do_user_streams = c_loc(ts_do_user_streams_ptr)
  inquire(iolength=transfer_struct_c%ts_do_user_streams_f_byte_size) fortran_type_f%ts_do_user_streams
#ifdef ifort
  transfer_struct_c%ts_do_user_streams_f_byte_size = transfer_struct_c%ts_do_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_additional_mvout = .FALSE.
  ts_do_additional_mvout_ptr => fortran_type_f%ts_do_additional_mvout
  transfer_struct_c%ts_do_additional_mvout = c_loc(ts_do_additional_mvout_ptr)
  inquire(iolength=transfer_struct_c%ts_do_additional_mvout_f_byte_size) fortran_type_f%ts_do_additional_mvout
#ifdef ifort
  transfer_struct_c%ts_do_additional_mvout_f_byte_size = transfer_struct_c%ts_do_additional_mvout_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_mvout_only = .FALSE.
  ts_do_mvout_only_ptr => fortran_type_f%ts_do_mvout_only
  transfer_struct_c%ts_do_mvout_only = c_loc(ts_do_mvout_only_ptr)
  inquire(iolength=transfer_struct_c%ts_do_mvout_only_f_byte_size) fortran_type_f%ts_do_mvout_only
#ifdef ifort
  transfer_struct_c%ts_do_mvout_only_f_byte_size = transfer_struct_c%ts_do_mvout_only_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_thermal_transonly = .FALSE.
  ts_do_thermal_transonly_ptr => fortran_type_f%ts_do_thermal_transonly
  transfer_struct_c%ts_do_thermal_transonly = c_loc(ts_do_thermal_transonly_ptr)
  inquire(iolength=transfer_struct_c%ts_do_thermal_transonly_f_byte_size) fortran_type_f%ts_do_thermal_transonly
#ifdef ifort
  transfer_struct_c%ts_do_thermal_transonly_f_byte_size = transfer_struct_c%ts_do_thermal_transonly_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_observation_geometry = .FALSE.
  ts_do_observation_geometry_ptr => fortran_type_f%ts_do_observation_geometry
  transfer_struct_c%ts_do_observation_geometry = c_loc(ts_do_observation_geometry_ptr)
  inquire(iolength=transfer_struct_c%ts_do_observation_geometry_f_byte_size) fortran_type_f%ts_do_observation_geometry
#ifdef ifort
  transfer_struct_c%ts_do_observation_geometry_f_byte_size = transfer_struct_c%ts_do_observation_geometry_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_doublet_geometry = .FALSE.
  ts_do_doublet_geometry_ptr => fortran_type_f%ts_do_doublet_geometry
  transfer_struct_c%ts_do_doublet_geometry = c_loc(ts_do_doublet_geometry_ptr)
  inquire(iolength=transfer_struct_c%ts_do_doublet_geometry_f_byte_size) fortran_type_f%ts_do_doublet_geometry
#ifdef ifort
  transfer_struct_c%ts_do_doublet_geometry_f_byte_size = transfer_struct_c%ts_do_doublet_geometry_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_boolean_c_init_only

subroutine lidort_modified_boolean_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_boolean

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_boolean_c_destroy

subroutine lidort_modified_boolean_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_boolean

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_boolean), pointer :: fortran_type_f_from
  type(lidort_modified_boolean), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_focorr = fortran_type_f_from%ts_do_focorr
  fortran_type_f_to%ts_do_focorr_external = fortran_type_f_from%ts_do_focorr_external
  fortran_type_f_to%ts_do_focorr_nadir = fortran_type_f_from%ts_do_focorr_nadir
  fortran_type_f_to%ts_do_focorr_outgoing = fortran_type_f_from%ts_do_focorr_outgoing
  fortran_type_f_to%ts_do_sscorr_truncation = fortran_type_f_from%ts_do_sscorr_truncation
  fortran_type_f_to%ts_do_sscorr_usephasfunc = fortran_type_f_from%ts_do_sscorr_usephasfunc
  fortran_type_f_to%ts_do_external_wleave = fortran_type_f_from%ts_do_external_wleave
  fortran_type_f_to%ts_do_double_convtest = fortran_type_f_from%ts_do_double_convtest
  fortran_type_f_to%ts_do_solar_sources = fortran_type_f_from%ts_do_solar_sources
  fortran_type_f_to%ts_do_refractive_geometry = fortran_type_f_from%ts_do_refractive_geometry
  fortran_type_f_to%ts_do_chapman_function = fortran_type_f_from%ts_do_chapman_function
  fortran_type_f_to%ts_do_rayleigh_only = fortran_type_f_from%ts_do_rayleigh_only
  fortran_type_f_to%ts_do_isotropic_only = fortran_type_f_from%ts_do_isotropic_only
  fortran_type_f_to%ts_do_no_azimuth = fortran_type_f_from%ts_do_no_azimuth
  fortran_type_f_to%ts_do_all_fourier = fortran_type_f_from%ts_do_all_fourier
  fortran_type_f_to%ts_do_deltam_scaling = fortran_type_f_from%ts_do_deltam_scaling
  fortran_type_f_to%ts_do_solution_saving = fortran_type_f_from%ts_do_solution_saving
  fortran_type_f_to%ts_do_bvp_telescoping = fortran_type_f_from%ts_do_bvp_telescoping
  fortran_type_f_to%ts_do_user_streams = fortran_type_f_from%ts_do_user_streams
  fortran_type_f_to%ts_do_additional_mvout = fortran_type_f_from%ts_do_additional_mvout
  fortran_type_f_to%ts_do_mvout_only = fortran_type_f_from%ts_do_mvout_only
  fortran_type_f_to%ts_do_thermal_transonly = fortran_type_f_from%ts_do_thermal_transonly
  fortran_type_f_to%ts_do_observation_geometry = fortran_type_f_from%ts_do_observation_geometry
  fortran_type_f_to%ts_do_doublet_geometry = fortran_type_f_from%ts_do_doublet_geometry
  

end subroutine lidort_modified_boolean_c_copy

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_control_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_control

  type(lidort_modified_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_control_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_control_c_alloc_init

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_control_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_nmoments_input_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_nmoments_input = 0
  ts_nmoments_input_ptr => fortran_type_f%ts_nmoments_input
  transfer_struct_c%ts_nmoments_input = c_loc(ts_nmoments_input_ptr)
  inquire(iolength=transfer_struct_c%ts_nmoments_input_f_byte_size) fortran_type_f%ts_nmoments_input
#ifdef ifort
  transfer_struct_c%ts_nmoments_input_f_byte_size = transfer_struct_c%ts_nmoments_input_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_control_c_init_only

subroutine lidort_modified_control_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_control

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_control_c_destroy

subroutine lidort_modified_control_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_control

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_control), pointer :: fortran_type_f_from
  type(lidort_modified_control), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_nmoments_input = fortran_type_f_from%ts_nmoments_input
  

end subroutine lidort_modified_control_c_copy

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_sunrays_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_sunrays

  type(lidort_modified_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_sunrays_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_sunrays_c_alloc_init

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_sunrays_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_nbeams_ptr
  real(c_double), dimension(:), pointer :: ts_beam_szas_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_nbeams = 0
  ts_nbeams_ptr => fortran_type_f%ts_nbeams
  transfer_struct_c%ts_nbeams = c_loc(ts_nbeams_ptr)
  inquire(iolength=transfer_struct_c%ts_nbeams_f_byte_size) fortran_type_f%ts_nbeams
#ifdef ifort
  transfer_struct_c%ts_nbeams_f_byte_size = transfer_struct_c%ts_nbeams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_beam_szas = 0_fpk
  ts_beam_szas_ptr => fortran_type_f%ts_beam_szas
  transfer_struct_c%ts_beam_szas = c_loc(ts_beam_szas_ptr(&
    lbound(fortran_type_f%ts_beam_szas,1)))
  inquire(iolength=transfer_struct_c%ts_beam_szas_f_byte_size) fortran_type_f%ts_beam_szas(&
    lbound(fortran_type_f%ts_beam_szas,1))
#ifdef ifort
  transfer_struct_c%ts_beam_szas_f_byte_size = transfer_struct_c%ts_beam_szas_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_beam_szas_f_shapes(1) = size(fortran_type_f%ts_beam_szas, 1)
  
  
end subroutine lidort_modified_sunrays_c_init_only

subroutine lidort_modified_sunrays_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_sunrays_c_destroy

subroutine lidort_modified_sunrays_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_sunrays), pointer :: fortran_type_f_from
  type(lidort_modified_sunrays), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_nbeams = fortran_type_f_from%ts_nbeams
  fortran_type_f_to%ts_beam_szas = fortran_type_f_from%ts_beam_szas
  

end subroutine lidort_modified_sunrays_c_copy

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_uservalues_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_uservalues

  type(lidort_modified_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_uservalues_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_uservalues_c_alloc_init

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_uservalues_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_n_user_relazms_ptr
  real(c_double), dimension(:), pointer :: ts_user_relazms_ptr
  integer(c_int), pointer :: ts_n_user_streams_ptr
  real(c_double), dimension(:), pointer :: ts_user_angles_input_ptr
  real(c_double), dimension(:), pointer :: ts_user_levels_ptr
  real(c_double), pointer :: ts_geometry_specheight_ptr
  integer(c_int), pointer :: ts_n_user_obsgeoms_ptr
  real(c_double), dimension(:,:), pointer :: ts_user_obsgeoms_input_ptr
  integer(c_int), pointer :: ts_n_user_doublets_ptr
  real(c_double), dimension(:,:), pointer :: ts_user_doublets_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_n_user_relazms = 0
  ts_n_user_relazms_ptr => fortran_type_f%ts_n_user_relazms
  transfer_struct_c%ts_n_user_relazms = c_loc(ts_n_user_relazms_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_relazms_f_byte_size) fortran_type_f%ts_n_user_relazms
#ifdef ifort
  transfer_struct_c%ts_n_user_relazms_f_byte_size = transfer_struct_c%ts_n_user_relazms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_user_relazms = 0_fpk
  ts_user_relazms_ptr => fortran_type_f%ts_user_relazms
  transfer_struct_c%ts_user_relazms = c_loc(ts_user_relazms_ptr(&
    lbound(fortran_type_f%ts_user_relazms,1)))
  inquire(iolength=transfer_struct_c%ts_user_relazms_f_byte_size) fortran_type_f%ts_user_relazms(&
    lbound(fortran_type_f%ts_user_relazms,1))
#ifdef ifort
  transfer_struct_c%ts_user_relazms_f_byte_size = transfer_struct_c%ts_user_relazms_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_relazms_f_shapes(1) = size(fortran_type_f%ts_user_relazms, 1)
  
  
  fortran_type_f%ts_n_user_streams = 0
  ts_n_user_streams_ptr => fortran_type_f%ts_n_user_streams
  transfer_struct_c%ts_n_user_streams = c_loc(ts_n_user_streams_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_streams_f_byte_size) fortran_type_f%ts_n_user_streams
#ifdef ifort
  transfer_struct_c%ts_n_user_streams_f_byte_size = transfer_struct_c%ts_n_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_user_angles_input = 0_fpk
  ts_user_angles_input_ptr => fortran_type_f%ts_user_angles_input
  transfer_struct_c%ts_user_angles_input = c_loc(ts_user_angles_input_ptr(&
    lbound(fortran_type_f%ts_user_angles_input,1)))
  inquire(iolength=transfer_struct_c%ts_user_angles_input_f_byte_size) fortran_type_f%ts_user_angles_input(&
    lbound(fortran_type_f%ts_user_angles_input,1))
#ifdef ifort
  transfer_struct_c%ts_user_angles_input_f_byte_size = transfer_struct_c%ts_user_angles_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_angles_input_f_shapes(1) = size(fortran_type_f%ts_user_angles_input, 1)
  
  
  fortran_type_f%ts_user_levels = 0_fpk
  ts_user_levels_ptr => fortran_type_f%ts_user_levels
  transfer_struct_c%ts_user_levels = c_loc(ts_user_levels_ptr(&
    lbound(fortran_type_f%ts_user_levels,1)))
  inquire(iolength=transfer_struct_c%ts_user_levels_f_byte_size) fortran_type_f%ts_user_levels(&
    lbound(fortran_type_f%ts_user_levels,1))
#ifdef ifort
  transfer_struct_c%ts_user_levels_f_byte_size = transfer_struct_c%ts_user_levels_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_levels_f_shapes(1) = size(fortran_type_f%ts_user_levels, 1)
  
  
  fortran_type_f%ts_geometry_specheight = 0_fpk
  ts_geometry_specheight_ptr => fortran_type_f%ts_geometry_specheight
  transfer_struct_c%ts_geometry_specheight = c_loc(ts_geometry_specheight_ptr)
  inquire(iolength=transfer_struct_c%ts_geometry_specheight_f_byte_size) fortran_type_f%ts_geometry_specheight
#ifdef ifort
  transfer_struct_c%ts_geometry_specheight_f_byte_size = transfer_struct_c%ts_geometry_specheight_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_n_user_obsgeoms = 0
  ts_n_user_obsgeoms_ptr => fortran_type_f%ts_n_user_obsgeoms
  transfer_struct_c%ts_n_user_obsgeoms = c_loc(ts_n_user_obsgeoms_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_obsgeoms_f_byte_size) fortran_type_f%ts_n_user_obsgeoms
#ifdef ifort
  transfer_struct_c%ts_n_user_obsgeoms_f_byte_size = transfer_struct_c%ts_n_user_obsgeoms_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_user_obsgeoms_input = 0_fpk
  ts_user_obsgeoms_input_ptr => fortran_type_f%ts_user_obsgeoms_input
  transfer_struct_c%ts_user_obsgeoms_input = c_loc(ts_user_obsgeoms_input_ptr(&
    lbound(fortran_type_f%ts_user_obsgeoms_input,1),&
    lbound(fortran_type_f%ts_user_obsgeoms_input,2)))
  inquire(iolength=transfer_struct_c%ts_user_obsgeoms_input_f_byte_size) fortran_type_f%ts_user_obsgeoms_input(&
    lbound(fortran_type_f%ts_user_obsgeoms_input,1),&
    lbound(fortran_type_f%ts_user_obsgeoms_input,2))
#ifdef ifort
  transfer_struct_c%ts_user_obsgeoms_input_f_byte_size = transfer_struct_c%ts_user_obsgeoms_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_obsgeoms_input_f_shapes(1) = size(fortran_type_f%ts_user_obsgeoms_input, 1)
  transfer_struct_c%ts_user_obsgeoms_input_f_shapes(2) = size(fortran_type_f%ts_user_obsgeoms_input, 2)
  
  
  fortran_type_f%ts_n_user_doublets = 0
  ts_n_user_doublets_ptr => fortran_type_f%ts_n_user_doublets
  transfer_struct_c%ts_n_user_doublets = c_loc(ts_n_user_doublets_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_doublets_f_byte_size) fortran_type_f%ts_n_user_doublets
#ifdef ifort
  transfer_struct_c%ts_n_user_doublets_f_byte_size = transfer_struct_c%ts_n_user_doublets_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_user_doublets = 0_fpk
  ts_user_doublets_ptr => fortran_type_f%ts_user_doublets
  transfer_struct_c%ts_user_doublets = c_loc(ts_user_doublets_ptr(&
    lbound(fortran_type_f%ts_user_doublets,1),&
    lbound(fortran_type_f%ts_user_doublets,2)))
  inquire(iolength=transfer_struct_c%ts_user_doublets_f_byte_size) fortran_type_f%ts_user_doublets(&
    lbound(fortran_type_f%ts_user_doublets,1),&
    lbound(fortran_type_f%ts_user_doublets,2))
#ifdef ifort
  transfer_struct_c%ts_user_doublets_f_byte_size = transfer_struct_c%ts_user_doublets_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_doublets_f_shapes(1) = size(fortran_type_f%ts_user_doublets, 1)
  transfer_struct_c%ts_user_doublets_f_shapes(2) = size(fortran_type_f%ts_user_doublets, 2)
  
  
end subroutine lidort_modified_uservalues_c_init_only

subroutine lidort_modified_uservalues_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_uservalues_c_destroy

subroutine lidort_modified_uservalues_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_uservalues), pointer :: fortran_type_f_from
  type(lidort_modified_uservalues), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_n_user_relazms = fortran_type_f_from%ts_n_user_relazms
  fortran_type_f_to%ts_user_relazms = fortran_type_f_from%ts_user_relazms
  fortran_type_f_to%ts_n_user_streams = fortran_type_f_from%ts_n_user_streams
  fortran_type_f_to%ts_user_angles_input = fortran_type_f_from%ts_user_angles_input
  fortran_type_f_to%ts_user_levels = fortran_type_f_from%ts_user_levels
  fortran_type_f_to%ts_geometry_specheight = fortran_type_f_from%ts_geometry_specheight
  fortran_type_f_to%ts_n_user_obsgeoms = fortran_type_f_from%ts_n_user_obsgeoms
  fortran_type_f_to%ts_user_obsgeoms_input = fortran_type_f_from%ts_user_obsgeoms_input
  fortran_type_f_to%ts_n_user_doublets = fortran_type_f_from%ts_n_user_doublets
  fortran_type_f_to%ts_user_doublets = fortran_type_f_from%ts_user_doublets
  

end subroutine lidort_modified_uservalues_c_copy

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_chapman_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_chapman

  type(lidort_modified_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_chapman_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_chapman_c_alloc_init

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_chapman_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  real(c_double), pointer :: ts_earth_radius_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_earth_radius = 0_fpk
  ts_earth_radius_ptr => fortran_type_f%ts_earth_radius
  transfer_struct_c%ts_earth_radius = c_loc(ts_earth_radius_ptr)
  inquire(iolength=transfer_struct_c%ts_earth_radius_f_byte_size) fortran_type_f%ts_earth_radius
#ifdef ifort
  transfer_struct_c%ts_earth_radius_f_byte_size = transfer_struct_c%ts_earth_radius_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_chapman_c_init_only

subroutine lidort_modified_chapman_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_chapman

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_chapman_c_destroy

subroutine lidort_modified_chapman_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_chapman

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_chapman), pointer :: fortran_type_f_from
  type(lidort_modified_chapman), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_earth_radius = fortran_type_f_from%ts_earth_radius
  

end subroutine lidort_modified_chapman_c_copy

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_optical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_optical

  type(lidort_modified_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_optical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_optical_c_alloc_init

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_optical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  real(c_double), dimension(:), pointer :: ts_omega_total_input_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_omega_total_input = 0_fpk
  ts_omega_total_input_ptr => fortran_type_f%ts_omega_total_input
  transfer_struct_c%ts_omega_total_input = c_loc(ts_omega_total_input_ptr(&
    lbound(fortran_type_f%ts_omega_total_input,1)))
  inquire(iolength=transfer_struct_c%ts_omega_total_input_f_byte_size) fortran_type_f%ts_omega_total_input(&
    lbound(fortran_type_f%ts_omega_total_input,1))
#ifdef ifort
  transfer_struct_c%ts_omega_total_input_f_byte_size = transfer_struct_c%ts_omega_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_omega_total_input_f_shapes(1) = size(fortran_type_f%ts_omega_total_input, 1)
  
  
end subroutine lidort_modified_optical_c_init_only

subroutine lidort_modified_optical_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_optical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_optical_c_destroy

subroutine lidort_modified_optical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_optical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_optical), pointer :: fortran_type_f_from
  type(lidort_modified_optical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_omega_total_input = fortran_type_f_from%ts_omega_total_input
  

end subroutine lidort_modified_optical_c_copy

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_inputs

  type(lidort_modified_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_inputs_c_alloc_init

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def_m
  use lidort_pars_m

  type(lidort_modified_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  type(lidort_modified_boolean), pointer :: mbool_ptr
  type(lidort_modified_control), pointer :: mcont_ptr
  type(lidort_modified_sunrays), pointer :: msunrays_ptr
  type(lidort_modified_uservalues), pointer :: muserval_ptr
  type(lidort_modified_chapman), pointer :: mchapman_ptr
  type(lidort_modified_optical), pointer :: moptical_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  mbool_ptr => fortran_type_f%mbool
  transfer_struct_c%mbool = c_loc(mbool_ptr)
  inquire(iolength=transfer_struct_c%mbool_f_byte_size) fortran_type_f%mbool
#ifdef ifort
  transfer_struct_c%mbool_f_byte_size = transfer_struct_c%mbool_f_byte_size * 4
#endif
  
  
  mcont_ptr => fortran_type_f%mcont
  transfer_struct_c%mcont = c_loc(mcont_ptr)
  inquire(iolength=transfer_struct_c%mcont_f_byte_size) fortran_type_f%mcont
#ifdef ifort
  transfer_struct_c%mcont_f_byte_size = transfer_struct_c%mcont_f_byte_size * 4
#endif
  
  
  msunrays_ptr => fortran_type_f%msunrays
  transfer_struct_c%msunrays = c_loc(msunrays_ptr)
  inquire(iolength=transfer_struct_c%msunrays_f_byte_size) fortran_type_f%msunrays
#ifdef ifort
  transfer_struct_c%msunrays_f_byte_size = transfer_struct_c%msunrays_f_byte_size * 4
#endif
  
  
  muserval_ptr => fortran_type_f%muserval
  transfer_struct_c%muserval = c_loc(muserval_ptr)
  inquire(iolength=transfer_struct_c%muserval_f_byte_size) fortran_type_f%muserval
#ifdef ifort
  transfer_struct_c%muserval_f_byte_size = transfer_struct_c%muserval_f_byte_size * 4
#endif
  
  
  mchapman_ptr => fortran_type_f%mchapman
  transfer_struct_c%mchapman = c_loc(mchapman_ptr)
  inquire(iolength=transfer_struct_c%mchapman_f_byte_size) fortran_type_f%mchapman
#ifdef ifort
  transfer_struct_c%mchapman_f_byte_size = transfer_struct_c%mchapman_f_byte_size * 4
#endif
  
  
  moptical_ptr => fortran_type_f%moptical
  transfer_struct_c%moptical = c_loc(moptical_ptr)
  inquire(iolength=transfer_struct_c%moptical_f_byte_size) fortran_type_f%moptical
#ifdef ifort
  transfer_struct_c%moptical_f_byte_size = transfer_struct_c%moptical_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_inputs_c_init_only

subroutine lidort_modified_inputs_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_inputs_c_destroy

subroutine lidort_modified_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def_m, only : lidort_modified_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_inputs), pointer :: fortran_type_f_from
  type(lidort_modified_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%mbool = fortran_type_f_from%mbool
  fortran_type_f_to%mcont = fortran_type_f_from%mcont
  fortran_type_f_to%msunrays = fortran_type_f_from%msunrays
  fortran_type_f_to%muserval = fortran_type_f_from%muserval
  fortran_type_f_to%mchapman = fortran_type_f_from%mchapman
  fortran_type_f_to%moptical = fortran_type_f_from%moptical
  

end subroutine lidort_modified_inputs_c_copy


! Wrapper for character variable "bs_brdf_names" of type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def_m" in file: "brdf_sup_inputs_def.F90"
subroutine brdf_sup_inputs_bs_brdf_names_get(fortran_type_c, bs_brdf_names_in_shape_1, &
      bs_brdf_names_in_len, &
      bs_brdf_names_in) bind(C)
  use brdf_sup_inputs_def_m, only : brdf_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: bs_brdf_names_in_shape_1
  integer(c_int), intent(in) :: bs_brdf_names_in_len
  character(kind=c_char) , intent(inout) :: bs_brdf_names_in(bs_brdf_names_in_shape_1, bs_brdf_names_in_len+1)

  type(brdf_sup_inputs), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%bs_brdf_names,1)
  do dim_idx_1 = 1, bs_brdf_names_in_shape_1
    do len_idx = 1, bs_brdf_names_in_len
      bs_brdf_names_in(dim_idx_1, len_idx) = &
          fortran_type_f%bs_brdf_names(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%bs_brdf_names(dim_idx_1-1+lb_1)(:))+1
    bs_brdf_names_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine brdf_sup_inputs_bs_brdf_names_get
! Wrapper for character variable "bs_inputmessages" of type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
subroutine brdf_input_exception_handling_bs_inputmessages_get(fortran_type_c, bs_inputmessages_in_shape_1, &
      bs_inputmessages_in_len, &
      bs_inputmessages_in) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: bs_inputmessages_in_shape_1
  integer(c_int), intent(in) :: bs_inputmessages_in_len
  character(kind=c_char) , intent(inout) :: bs_inputmessages_in(bs_inputmessages_in_shape_1, bs_inputmessages_in_len+1)

  type(brdf_input_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%bs_inputmessages,1)
  do dim_idx_1 = 1, bs_inputmessages_in_shape_1
    do len_idx = 1, bs_inputmessages_in_len
      bs_inputmessages_in(dim_idx_1, len_idx) = &
          fortran_type_f%bs_inputmessages(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%bs_inputmessages(dim_idx_1-1+lb_1)(:))+1
    bs_inputmessages_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine brdf_input_exception_handling_bs_inputmessages_get
! Wrapper for character variable "bs_inputactions" of type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
subroutine brdf_input_exception_handling_bs_inputactions_get(fortran_type_c, bs_inputactions_in_shape_1, &
      bs_inputactions_in_len, &
      bs_inputactions_in) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: bs_inputactions_in_shape_1
  integer(c_int), intent(in) :: bs_inputactions_in_len
  character(kind=c_char) , intent(inout) :: bs_inputactions_in(bs_inputactions_in_shape_1, bs_inputactions_in_len+1)

  type(brdf_input_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%bs_inputactions,1)
  do dim_idx_1 = 1, bs_inputactions_in_shape_1
    do len_idx = 1, bs_inputactions_in_len
      bs_inputactions_in(dim_idx_1, len_idx) = &
          fortran_type_f%bs_inputactions(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%bs_inputactions(dim_idx_1-1+lb_1)(:))+1
    bs_inputactions_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine brdf_input_exception_handling_bs_inputactions_get
! Wrapper for character variable "bs_outputmessages" of type: "brdf_output_exception_handling" from module: "brdf_sup_outputs_def_m" in file: "brdf_sup_outputs_def.F90"
subroutine brdf_output_exception_handling_bs_outputmessages_get(fortran_type_c, bs_outputmessages_in_shape_1, &
      bs_outputmessages_in_len, &
      bs_outputmessages_in) bind(C)
  use brdf_sup_outputs_def_m, only : brdf_output_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: bs_outputmessages_in_shape_1
  integer(c_int), intent(in) :: bs_outputmessages_in_len
  character(kind=c_char) , intent(inout) :: bs_outputmessages_in(bs_outputmessages_in_shape_1, bs_outputmessages_in_len+1)

  type(brdf_output_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%bs_outputmessages,1)
  do dim_idx_1 = 1, bs_outputmessages_in_shape_1
    do len_idx = 1, bs_outputmessages_in_len
      bs_outputmessages_in(dim_idx_1, len_idx) = &
          fortran_type_f%bs_outputmessages(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%bs_outputmessages(dim_idx_1-1+lb_1)(:))+1
    bs_outputmessages_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine brdf_output_exception_handling_bs_outputmessages_get
! Wrapper for character variable "sl_sleave_datapath" of type: "sleave_sup_inputs" from module: "sleave_sup_inputs_def_m" in file: "sleave_sup_inputs_def.F90"
subroutine sleave_sup_inputs_sl_sleave_datapath_get(fortran_type_c, sl_sleave_datapath_in_len, &
      sl_sleave_datapath_in) bind(C)
  use sleave_sup_inputs_def_m, only : sleave_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: sl_sleave_datapath_in_len
  character(kind=c_char) , intent(inout) :: sl_sleave_datapath_in(sl_sleave_datapath_in_len+1)

  type(sleave_sup_inputs), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, sl_sleave_datapath_in_len
    sl_sleave_datapath_in(len_idx) = &
      fortran_type_f%sl_sleave_datapath(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%sl_sleave_datapath(:))+1
  sl_sleave_datapath_in(len_idx) = c_null_char

end subroutine sleave_sup_inputs_sl_sleave_datapath_get
! Wrapper for character variable "ts_columnwf_names" of type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
subroutine fixed_lincontrol_ts_columnwf_names_get(fortran_type_c, ts_columnwf_names_in_shape_1, &
      ts_columnwf_names_in_len, &
      ts_columnwf_names_in) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_columnwf_names_in_shape_1
  integer(c_int), intent(in) :: ts_columnwf_names_in_len
  character(kind=c_char) , intent(inout) :: ts_columnwf_names_in(ts_columnwf_names_in_shape_1, ts_columnwf_names_in_len+1)

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_columnwf_names,1)
  do dim_idx_1 = 1, ts_columnwf_names_in_shape_1
    do len_idx = 1, ts_columnwf_names_in_len
      ts_columnwf_names_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_columnwf_names(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_columnwf_names(dim_idx_1-1+lb_1)(:))+1
    ts_columnwf_names_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine fixed_lincontrol_ts_columnwf_names_get
! Wrapper for character variable "ts_profilewf_names" of type: "lidort_fixed_lincontrol" from module: "lidort_lin_inputs_def_m" in file: "lidort_lin_inputs_def.F90"
subroutine fixed_lincontrol_ts_profilewf_names_get(fortran_type_c, ts_profilewf_names_in_shape_1, &
      ts_profilewf_names_in_len, &
      ts_profilewf_names_in) bind(C)
  use lidort_lin_inputs_def_m, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_profilewf_names_in_shape_1
  integer(c_int), intent(in) :: ts_profilewf_names_in_len
  character(kind=c_char) , intent(inout) :: ts_profilewf_names_in(ts_profilewf_names_in_shape_1, ts_profilewf_names_in_len+1)

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_profilewf_names,1)
  do dim_idx_1 = 1, ts_profilewf_names_in_shape_1
    do len_idx = 1, ts_profilewf_names_in_len
      ts_profilewf_names_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_profilewf_names(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_profilewf_names(dim_idx_1-1+lb_1)(:))+1
    ts_profilewf_names_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine fixed_lincontrol_ts_profilewf_names_get
! Wrapper for character variable "ts_checkmessages" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_checkmessages_get(fortran_type_c, ts_checkmessages_in_shape_1, &
      ts_checkmessages_in_len, &
      ts_checkmessages_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_checkmessages_in_shape_1
  integer(c_int), intent(in) :: ts_checkmessages_in_len
  character(kind=c_char) , intent(inout) :: ts_checkmessages_in(ts_checkmessages_in_shape_1, ts_checkmessages_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_checkmessages,1)
  do dim_idx_1 = 1, ts_checkmessages_in_shape_1
    do len_idx = 1, ts_checkmessages_in_len
      ts_checkmessages_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_checkmessages(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_checkmessages(dim_idx_1-1+lb_1)(:))+1
    ts_checkmessages_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine exception_handling_ts_checkmessages_get
! Wrapper for character variable "ts_actions" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_actions_get(fortran_type_c, ts_actions_in_shape_1, &
      ts_actions_in_len, &
      ts_actions_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_actions_in_shape_1
  integer(c_int), intent(in) :: ts_actions_in_len
  character(kind=c_char) , intent(inout) :: ts_actions_in(ts_actions_in_shape_1, ts_actions_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_actions,1)
  do dim_idx_1 = 1, ts_actions_in_shape_1
    do len_idx = 1, ts_actions_in_len
      ts_actions_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_actions(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_actions(dim_idx_1-1+lb_1)(:))+1
    ts_actions_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine exception_handling_ts_actions_get
! Wrapper for character variable "ts_message" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_message_get(fortran_type_c, ts_message_in_len, &
      ts_message_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_message_in_len
  character(kind=c_char) , intent(inout) :: ts_message_in(ts_message_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_message_in_len
    ts_message_in(len_idx) = &
      fortran_type_f%ts_message(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_message(:))+1
  ts_message_in(len_idx) = c_null_char

end subroutine exception_handling_ts_message_get
! Wrapper for character variable "ts_trace_1" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_1_get(fortran_type_c, ts_trace_1_in_len, &
      ts_trace_1_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_trace_1_in_len
  character(kind=c_char) , intent(inout) :: ts_trace_1_in(ts_trace_1_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_trace_1_in_len
    ts_trace_1_in(len_idx) = &
      fortran_type_f%ts_trace_1(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_trace_1(:))+1
  ts_trace_1_in(len_idx) = c_null_char

end subroutine exception_handling_ts_trace_1_get
! Wrapper for character variable "ts_trace_2" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_2_get(fortran_type_c, ts_trace_2_in_len, &
      ts_trace_2_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_trace_2_in_len
  character(kind=c_char) , intent(inout) :: ts_trace_2_in(ts_trace_2_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_trace_2_in_len
    ts_trace_2_in(len_idx) = &
      fortran_type_f%ts_trace_2(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_trace_2(:))+1
  ts_trace_2_in(len_idx) = c_null_char

end subroutine exception_handling_ts_trace_2_get
! Wrapper for character variable "ts_trace_3" of type: "lidort_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_3_get(fortran_type_c, ts_trace_3_in_len, &
      ts_trace_3_in) bind(C)
  use lidort_outputs_def_m, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_trace_3_in_len
  character(kind=c_char) , intent(inout) :: ts_trace_3_in(ts_trace_3_in_len+1)

  type(lidort_exception_handling), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_trace_3_in_len
    ts_trace_3_in(len_idx) = &
      fortran_type_f%ts_trace_3(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_trace_3(:))+1
  ts_trace_3_in(len_idx) = c_null_char

end subroutine exception_handling_ts_trace_3_get
! Wrapper for character variable "ts_inputmessages" of type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine input_exception_handling_ts_inputmessages_get(fortran_type_c, ts_inputmessages_in_shape_1, &
      ts_inputmessages_in_len, &
      ts_inputmessages_in) bind(C)
  use lidort_outputs_def_m, only : lidort_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_inputmessages_in_shape_1
  integer(c_int), intent(in) :: ts_inputmessages_in_len
  character(kind=c_char) , intent(inout) :: ts_inputmessages_in(ts_inputmessages_in_shape_1, ts_inputmessages_in_len+1)

  type(lidort_input_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_inputmessages,1)
  do dim_idx_1 = 1, ts_inputmessages_in_shape_1
    do len_idx = 1, ts_inputmessages_in_len
      ts_inputmessages_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_inputmessages(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_inputmessages(dim_idx_1-1+lb_1)(:))+1
    ts_inputmessages_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine input_exception_handling_ts_inputmessages_get
! Wrapper for character variable "ts_inputactions" of type: "lidort_input_exception_handling" from module: "lidort_outputs_def_m" in file: "lidort_outputs_def.F90"
subroutine input_exception_handling_ts_inputactions_get(fortran_type_c, ts_inputactions_in_shape_1, &
      ts_inputactions_in_len, &
      ts_inputactions_in) bind(C)
  use lidort_outputs_def_m, only : lidort_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_inputactions_in_shape_1
  integer(c_int), intent(in) :: ts_inputactions_in_len
  character(kind=c_char) , intent(inout) :: ts_inputactions_in(ts_inputactions_in_shape_1, ts_inputactions_in_len+1)

  type(lidort_input_exception_handling), pointer :: fortran_type_f
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  lb_1 = lbound(fortran_type_f%ts_inputactions,1)
  do dim_idx_1 = 1, ts_inputactions_in_shape_1
    do len_idx = 1, ts_inputactions_in_len
      ts_inputactions_in(dim_idx_1, len_idx) = &
          fortran_type_f%ts_inputactions(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(fortran_type_f%ts_inputactions(dim_idx_1-1+lb_1)(:))+1
    ts_inputactions_in(dim_idx_1, len_idx) = c_null_char
  end do

end subroutine input_exception_handling_ts_inputactions_get
! Wrapper for character variable "ts_input_write_filename" of type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
subroutine fixed_write_ts_input_write_filename_get(fortran_type_c, ts_input_write_filename_in_len, &
      ts_input_write_filename_in) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_input_write_filename_in_len
  character(kind=c_char) , intent(inout) :: ts_input_write_filename_in(ts_input_write_filename_in_len+1)

  type(lidort_fixed_write), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_input_write_filename_in_len
    ts_input_write_filename_in(len_idx) = &
      fortran_type_f%ts_input_write_filename(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_input_write_filename(:))+1
  ts_input_write_filename_in(len_idx) = c_null_char

end subroutine fixed_write_ts_input_write_filename_get
! Wrapper for character variable "ts_scenario_write_filename" of type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
subroutine fixed_write_ts_scenario_write_filename_get(fortran_type_c, ts_scenario_write_filename_in_len, &
      ts_scenario_write_filename_in) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_scenario_write_filename_in_len
  character(kind=c_char) , intent(inout) :: ts_scenario_write_filename_in(ts_scenario_write_filename_in_len+1)

  type(lidort_fixed_write), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_scenario_write_filename_in_len
    ts_scenario_write_filename_in(len_idx) = &
      fortran_type_f%ts_scenario_write_filename(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_scenario_write_filename(:))+1
  ts_scenario_write_filename_in(len_idx) = c_null_char

end subroutine fixed_write_ts_scenario_write_filename_get
! Wrapper for character variable "ts_fourier_write_filename" of type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
subroutine fixed_write_ts_fourier_write_filename_get(fortran_type_c, ts_fourier_write_filename_in_len, &
      ts_fourier_write_filename_in) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_fourier_write_filename_in_len
  character(kind=c_char) , intent(inout) :: ts_fourier_write_filename_in(ts_fourier_write_filename_in_len+1)

  type(lidort_fixed_write), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_fourier_write_filename_in_len
    ts_fourier_write_filename_in(len_idx) = &
      fortran_type_f%ts_fourier_write_filename(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_fourier_write_filename(:))+1
  ts_fourier_write_filename_in(len_idx) = c_null_char

end subroutine fixed_write_ts_fourier_write_filename_get
! Wrapper for character variable "ts_results_write_filename" of type: "lidort_fixed_write" from module: "lidort_inputs_def_m" in file: "lidort_inputs_def.F90"
subroutine fixed_write_ts_results_write_filename_get(fortran_type_c, ts_results_write_filename_in_len, &
      ts_results_write_filename_in) bind(C)
  use lidort_inputs_def_m, only : lidort_fixed_write

  type(c_ptr), intent(inout) :: fortran_type_c
  integer(c_int), intent(in) :: ts_results_write_filename_in_len
  character(kind=c_char) , intent(inout) :: ts_results_write_filename_in(ts_results_write_filename_in_len+1)

  type(lidort_fixed_write), pointer :: fortran_type_f
  integer :: len_idx

  call c_f_pointer(fortran_type_c, fortran_type_f)

  do len_idx = 1, ts_results_write_filename_in_len
    ts_results_write_filename_in(len_idx) = &
      fortran_type_f%ts_results_write_filename(len_idx:len_idx)
  end do
  len_idx = len_trim(fortran_type_f%ts_results_write_filename(:))+1
  ts_results_write_filename_in(len_idx) = c_null_char

end subroutine fixed_write_ts_results_write_filename_get


end module lidort_interface_types