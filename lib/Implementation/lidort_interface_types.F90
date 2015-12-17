module lidort_interface_types

use iso_c_binding
use lidort_pars

! This module was auto-generated 

implicit none


! Links to module: "lidort_pars" in file: "lidort_pars.F90.in"
type, bind(c) :: lidort_pars_c
  character(kind=c_char, len=3) :: lidort_version_number
  integer(c_int) :: lidort_inunit
  integer(c_int) :: lidort_scenunit
  integer(c_int) :: lidort_funit
  integer(c_int) :: lidort_resunit
  integer(c_int) :: lidort_errunit
  integer(c_int) :: lidort_dbgunit
  integer(c_int) :: max_messages
  integer(c_int) :: maxthreads
  integer(c_int) :: maxstreams
  integer(c_int) :: maxlayers
  integer(c_int) :: maxfinelayers
  integer(c_int) :: maxmoments_input
  integer(c_int) :: max_thermal_coeffs
  integer(c_int) :: maxbeams
  integer(c_int) :: max_user_streams
  integer(c_int) :: max_user_relazms
  integer(c_int) :: max_user_levels
  integer(c_int) :: max_partlayers
  integer(c_int) :: max_directions
  integer(c_int) :: max_brdf_kernels
  integer(c_int) :: max_brdf_parameters
  integer(c_int) :: maxstreams_brdf
  integer(c_int) :: max_msrs_muquad
  integer(c_int) :: max_msrs_phiquad
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
  integer(c_int) :: breonveg_idx
  integer(c_int) :: breonsoil_idx
  integer(c_int) :: maxbrdf_idx
  
end type lidort_pars_c



! Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
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

  
end type brdf_linsup_inputs_c

! Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
type, bind(c) :: brdf_linsup_outputs_c
  type(c_ptr) :: bs_ls_exactdb_brdfunc ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(4) :: bs_ls_exactdb_brdfunc_f_shapes
  integer(c_int) :: bs_ls_exactdb_brdfunc_f_byte_size

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

  type(c_ptr) :: bs_ls_emissivity ! dimension(MAX_SURFACEWFS, MAXSTREAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: bs_ls_emissivity_f_shapes
  integer(c_int) :: bs_ls_emissivity_f_byte_size

  type(c_ptr) :: bs_ls_user_emissivity ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: bs_ls_user_emissivity_f_shapes
  integer(c_int) :: bs_ls_user_emissivity_f_byte_size

  
end type brdf_linsup_outputs_c

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
type, bind(c) :: brdf_sup_inputs_c
  type(c_ptr) :: bs_do_user_streams ! scalar
  integer(c_int) :: bs_do_user_streams_f_byte_size

  type(c_ptr) :: bs_do_brdf_surface ! scalar
  integer(c_int) :: bs_do_brdf_surface_f_byte_size

  type(c_ptr) :: bs_do_surface_emission ! scalar
  integer(c_int) :: bs_do_surface_emission_f_byte_size

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

  type(c_ptr) :: bs_do_exactonly ! scalar
  integer(c_int) :: bs_do_exactonly_f_byte_size

  type(c_ptr) :: bs_do_glitter_msrcorr ! scalar
  integer(c_int) :: bs_do_glitter_msrcorr_f_byte_size

  type(c_ptr) :: bs_do_glitter_msrcorr_exactonly ! scalar
  integer(c_int) :: bs_do_glitter_msrcorr_exactonly_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_order ! scalar
  integer(c_int) :: bs_glitter_msrcorr_order_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_nmuquad ! scalar
  integer(c_int) :: bs_glitter_msrcorr_nmuquad_f_byte_size

  type(c_ptr) :: bs_glitter_msrcorr_nphiquad ! scalar
  integer(c_int) :: bs_glitter_msrcorr_nphiquad_f_byte_size

  
end type brdf_sup_inputs_c

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
type, bind(c) :: brdf_sup_outputs_c
  type(c_ptr) :: bs_exactdb_brdfunc ! dimension(MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS)
  integer(c_int), dimension(3) :: bs_exactdb_brdfunc_f_shapes
  integer(c_int) :: bs_exactdb_brdfunc_f_byte_size

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

  type(c_ptr) :: bs_emissivity ! dimension(MAXSTREAMS, MAXTHREADS)
  integer(c_int), dimension(2) :: bs_emissivity_f_shapes
  integer(c_int) :: bs_emissivity_f_byte_size

  type(c_ptr) :: bs_user_emissivity ! dimension(MAX_USER_STREAMS, MAXTHREADS)
  integer(c_int), dimension(2) :: bs_user_emissivity_f_shapes
  integer(c_int) :: bs_user_emissivity_f_byte_size

  
end type brdf_sup_outputs_c

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
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

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_lincontrol_c
  type(c_ptr) :: ts_do_column_linearization ! scalar
  integer(c_int) :: ts_do_column_linearization_f_byte_size

  type(c_ptr) :: ts_do_profile_linearization ! scalar
  integer(c_int) :: ts_do_profile_linearization_f_byte_size

  type(c_ptr) :: ts_do_surface_linearization ! scalar
  integer(c_int) :: ts_do_surface_linearization_f_byte_size

  type(c_ptr) :: ts_do_sleave_wfs ! scalar
  integer(c_int) :: ts_do_sleave_wfs_f_byte_size

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

  
end type lidort_fixed_lincontrol_c

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_linoptical_c
  type(c_ptr) :: ts_l_deltau_vert_input ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_l_deltau_vert_input_f_shapes
  integer(c_int) :: ts_l_deltau_vert_input_f_byte_size

  type(c_ptr) :: ts_l_omega_total_input ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_l_omega_total_input_f_shapes
  integer(c_int) :: ts_l_omega_total_input_f_byte_size

  type(c_ptr) :: ts_l_phasmoms_total_input ! dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(4) :: ts_l_phasmoms_total_input_f_shapes
  integer(c_int) :: ts_l_phasmoms_total_input_f_byte_size

  
end type lidort_fixed_linoptical_c

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_fixed_lininputs_c
  type(c_ptr) :: cont ! scalar
  integer(c_int) :: cont_f_byte_size

  type(c_ptr) :: optical ! scalar
  integer(c_int) :: optical_f_byte_size

  
end type lidort_fixed_lininputs_c

! Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
type, bind(c) :: lidort_modified_lininputs_c
  type(c_ptr) :: dummy ! scalar
  integer(c_int) :: dummy_f_byte_size

  
end type lidort_modified_lininputs_c

! Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linatmos_c
  type(c_ptr) :: ts_columnwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_columnwf_f_shapes
  integer(c_int) :: ts_columnwf_f_byte_size

  type(c_ptr) :: ts_mint_columnwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_mint_columnwf_f_shapes
  integer(c_int) :: ts_mint_columnwf_f_byte_size

  type(c_ptr) :: ts_flux_columnwf ! dimension(MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_flux_columnwf_f_shapes
  integer(c_int) :: ts_flux_columnwf_f_byte_size

  type(c_ptr) :: ts_profilewf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(6) :: ts_profilewf_f_shapes
  integer(c_int) :: ts_profilewf_f_byte_size

  type(c_ptr) :: ts_mint_profilewf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(6) :: ts_mint_profilewf_f_shapes
  integer(c_int) :: ts_mint_profilewf_f_byte_size

  type(c_ptr) :: ts_flux_profilewf ! dimension(MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(6) :: ts_flux_profilewf_f_shapes
  integer(c_int) :: ts_flux_profilewf_f_byte_size

  
end type lidort_linatmos_c

! Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linsurf_c
  type(c_ptr) :: ts_surfacewf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_surfacewf_f_shapes
  integer(c_int) :: ts_surfacewf_f_byte_size

  type(c_ptr) :: ts_mint_surfacewf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_mint_surfacewf_f_shapes
  integer(c_int) :: ts_mint_surfacewf_f_byte_size

  type(c_ptr) :: ts_flux_surfacewf ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(5) :: ts_flux_surfacewf_f_shapes
  integer(c_int) :: ts_flux_surfacewf_f_byte_size

  
end type lidort_linsurf_c

! Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
type, bind(c) :: lidort_linoutputs_c
  type(c_ptr) :: atmos ! scalar
  integer(c_int) :: atmos_f_byte_size

  type(c_ptr) :: surf ! scalar
  integer(c_int) :: surf_f_byte_size

  
end type lidort_linoutputs_c

! Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
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

  type(c_ptr) :: ts_ls_emissivity ! dimension(MAX_SURFACEWFS, MAXSTREAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_ls_emissivity_f_shapes
  integer(c_int) :: ts_ls_emissivity_f_byte_size

  type(c_ptr) :: ts_ls_user_emissivity ! dimension(MAX_SURFACEWFS, MAX_USER_STREAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_ls_user_emissivity_f_shapes
  integer(c_int) :: ts_ls_user_emissivity_f_byte_size

  
end type lidort_linsup_brdf_c

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
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

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
type, bind(c) :: lidort_linsup_ss_surf_c
  type(c_ptr) :: ts_surfacewf_db ! dimension(MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(3) :: ts_surfacewf_db_f_shapes
  integer(c_int) :: ts_surfacewf_db_f_byte_size

  
end type lidort_linsup_ss_surf_c

! Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
type, bind(c) :: lidort_linsup_ss_c
  type(c_ptr) :: atmos ! scalar
  integer(c_int) :: atmos_f_byte_size

  type(c_ptr) :: surf ! scalar
  integer(c_int) :: surf_f_byte_size

  
end type lidort_linsup_ss_c

! Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
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

! Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
type, bind(c) :: lidort_linsup_inout_c
  type(c_ptr) :: brdf ! scalar
  integer(c_int) :: brdf_f_byte_size

  type(c_ptr) :: ss ! scalar
  integer(c_int) :: ss_f_byte_size

  type(c_ptr) :: sleave ! scalar
  integer(c_int) :: sleave_f_byte_size

  
end type lidort_linsup_inout_c

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_main_outputs_c
  type(c_ptr) :: ts_intensity ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(4) :: ts_intensity_f_shapes
  integer(c_int) :: ts_intensity_f_byte_size

  type(c_ptr) :: ts_mean_intensity ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(4) :: ts_mean_intensity_f_shapes
  integer(c_int) :: ts_mean_intensity_f_byte_size

  type(c_ptr) :: ts_flux_integral ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS)
  integer(c_int), dimension(4) :: ts_flux_integral_f_shapes
  integer(c_int) :: ts_flux_integral_f_byte_size

  type(c_ptr) :: ts_dnflux_direct ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_dnflux_direct_f_shapes
  integer(c_int) :: ts_dnflux_direct_f_byte_size

  type(c_ptr) :: ts_dnmean_direct ! dimension(MAX_USER_LEVELS, MAXBEAMS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_dnmean_direct_f_shapes
  integer(c_int) :: ts_dnmean_direct_f_byte_size

  type(c_ptr) :: ts_fourier_saved ! dimension(MAXBEAMS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_fourier_saved_f_shapes
  integer(c_int) :: ts_fourier_saved_f_byte_size

  type(c_ptr) :: ts_n_geometries ! scalar
  integer(c_int) :: ts_n_geometries_f_byte_size

  
end type lidort_main_outputs_c

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
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

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
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

! Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
type, bind(c) :: lidort_outputs_c
  type(c_ptr) :: main ! scalar
  integer(c_int) :: main_f_byte_size

  type(c_ptr) :: status ! scalar
  integer(c_int) :: status_f_byte_size

  
end type lidort_outputs_c

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
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

  type(c_ptr) :: ts_emissivity ! dimension(MAXSTREAMS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_emissivity_f_shapes
  integer(c_int) :: ts_emissivity_f_byte_size

  type(c_ptr) :: ts_user_emissivity ! dimension(MAX_USER_STREAMS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_user_emissivity_f_shapes
  integer(c_int) :: ts_user_emissivity_f_byte_size

  
end type lidort_sup_brdf_c

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
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

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
type, bind(c) :: lidort_sup_ss_c
  type(c_ptr) :: ts_intensity_ss ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)
  integer(c_int), dimension(3) :: ts_intensity_ss_f_shapes
  integer(c_int) :: ts_intensity_ss_f_byte_size

  type(c_ptr) :: ts_intensity_db ! dimension(MAX_USER_LEVELS, MAX_GEOMETRIES)
  integer(c_int), dimension(2) :: ts_intensity_db_f_shapes
  integer(c_int) :: ts_intensity_db_f_byte_size

  
end type lidort_sup_ss_c

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
type, bind(c) :: lidort_sup_inout_c
  type(c_ptr) :: brdf ! scalar
  integer(c_int) :: brdf_f_byte_size

  type(c_ptr) :: ss ! scalar
  integer(c_int) :: ss_f_byte_size

  type(c_ptr) :: sleave ! scalar
  integer(c_int) :: sleave_f_byte_size

  
end type lidort_sup_inout_c

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_boolean_c
  type(c_ptr) :: ts_do_fullrad_mode ! scalar
  integer(c_int) :: ts_do_fullrad_mode_f_byte_size

  type(c_ptr) :: ts_do_sscorr_truncation ! scalar
  integer(c_int) :: ts_do_sscorr_truncation_f_byte_size

  type(c_ptr) :: ts_do_ss_external ! scalar
  integer(c_int) :: ts_do_ss_external_f_byte_size

  type(c_ptr) :: ts_do_ssfull ! scalar
  integer(c_int) :: ts_do_ssfull_f_byte_size

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

  type(c_ptr) :: ts_do_surface_leaving ! scalar
  integer(c_int) :: ts_do_surface_leaving_f_byte_size

  type(c_ptr) :: ts_do_sl_isotropic ! scalar
  integer(c_int) :: ts_do_sl_isotropic_f_byte_size

  
end type lidort_fixed_boolean_c

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_control_c
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

  
end type lidort_fixed_control_c

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_sunrays_c
  type(c_ptr) :: ts_flux_factor ! scalar
  integer(c_int) :: ts_flux_factor_f_byte_size

  
end type lidort_fixed_sunrays_c

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_uservalues_c
  type(c_ptr) :: ts_n_user_streams ! scalar
  integer(c_int) :: ts_n_user_streams_f_byte_size

  type(c_ptr) :: ts_n_user_levels ! scalar
  integer(c_int) :: ts_n_user_levels_f_byte_size

  
end type lidort_fixed_uservalues_c

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
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

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_fixed_optical_c
  type(c_ptr) :: ts_deltau_vert_input ! dimension(MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_deltau_vert_input_f_shapes
  integer(c_int) :: ts_deltau_vert_input_f_byte_size

  type(c_ptr) :: ts_phasmoms_total_input ! dimension(0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(3) :: ts_phasmoms_total_input_f_shapes
  integer(c_int) :: ts_phasmoms_total_input_f_byte_size

  type(c_ptr) :: ts_thermal_bb_input ! dimension(0:MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_thermal_bb_input_f_shapes
  integer(c_int) :: ts_thermal_bb_input_f_byte_size

  type(c_ptr) :: ts_lambertian_albedo ! dimension(MAXTHREADS)
  integer(c_int), dimension(1) :: ts_lambertian_albedo_f_shapes
  integer(c_int) :: ts_lambertian_albedo_f_byte_size

  type(c_ptr) :: ts_surface_bb_input ! dimension(MAXTHREADS)
  integer(c_int), dimension(1) :: ts_surface_bb_input_f_shapes
  integer(c_int) :: ts_surface_bb_input_f_byte_size

  
end type lidort_fixed_optical_c

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
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

  
end type lidort_fixed_inputs_c

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_boolean_c
  type(c_ptr) :: ts_do_sscorr_nadir ! scalar
  integer(c_int) :: ts_do_sscorr_nadir_f_byte_size

  type(c_ptr) :: ts_do_sscorr_outgoing ! scalar
  integer(c_int) :: ts_do_sscorr_outgoing_f_byte_size

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

  
end type lidort_modified_boolean_c

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_control_c
  type(c_ptr) :: ts_nmoments_input ! scalar
  integer(c_int) :: ts_nmoments_input_f_byte_size

  
end type lidort_modified_control_c

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_sunrays_c
  type(c_ptr) :: ts_nbeams ! scalar
  integer(c_int) :: ts_nbeams_f_byte_size

  type(c_ptr) :: ts_beam_szas ! dimension(MAXBEAMS)
  integer(c_int), dimension(1) :: ts_beam_szas_f_shapes
  integer(c_int) :: ts_beam_szas_f_byte_size

  
end type lidort_modified_sunrays_c

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_uservalues_c
  type(c_ptr) :: ts_n_user_relazms ! scalar
  integer(c_int) :: ts_n_user_relazms_f_byte_size

  type(c_ptr) :: ts_user_relazms ! dimension(MAX_USER_RELAZMS)
  integer(c_int), dimension(1) :: ts_user_relazms_f_shapes
  integer(c_int) :: ts_user_relazms_f_byte_size

  type(c_ptr) :: ts_user_angles_input ! dimension(MAX_USER_STREAMS)
  integer(c_int), dimension(1) :: ts_user_angles_input_f_shapes
  integer(c_int) :: ts_user_angles_input_f_byte_size

  type(c_ptr) :: ts_user_levels ! dimension(MAX_USER_LEVELS)
  integer(c_int), dimension(1) :: ts_user_levels_f_shapes
  integer(c_int) :: ts_user_levels_f_byte_size

  type(c_ptr) :: ts_geometry_specheight ! scalar
  integer(c_int) :: ts_geometry_specheight_f_byte_size

  
end type lidort_modified_uservalues_c

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_chapman_c
  type(c_ptr) :: ts_earth_radius ! scalar
  integer(c_int) :: ts_earth_radius_f_byte_size

  
end type lidort_modified_chapman_c

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
type, bind(c) :: lidort_modified_optical_c
  type(c_ptr) :: ts_omega_total_input ! dimension(MAXLAYERS, MAXTHREADS)
  integer(c_int), dimension(2) :: ts_omega_total_input_f_shapes
  integer(c_int) :: ts_omega_total_input_f_byte_size

  
end type lidort_modified_optical_c

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
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
  type(lidort_pars_c) :: pars_struct

  pars_struct%lidort_version_number = LIDORT_VERSION_NUMBER
  pars_struct%lidort_inunit = LIDORT_INUNIT
  pars_struct%lidort_scenunit = LIDORT_SCENUNIT
  pars_struct%lidort_funit = LIDORT_FUNIT
  pars_struct%lidort_resunit = LIDORT_RESUNIT
  pars_struct%lidort_errunit = LIDORT_ERRUNIT
  pars_struct%lidort_dbgunit = LIDORT_DBGUNIT
  pars_struct%max_messages = MAX_MESSAGES
  pars_struct%maxthreads = MAXTHREADS
  pars_struct%maxstreams = MAXSTREAMS
  pars_struct%maxlayers = MAXLAYERS
  pars_struct%maxfinelayers = MAXFINELAYERS
  pars_struct%maxmoments_input = MAXMOMENTS_INPUT
  pars_struct%max_thermal_coeffs = MAX_THERMAL_COEFFS
  pars_struct%maxbeams = MAXBEAMS
  pars_struct%max_user_streams = MAX_USER_STREAMS
  pars_struct%max_user_relazms = MAX_USER_RELAZMS
  pars_struct%max_user_levels = MAX_USER_LEVELS
  pars_struct%max_partlayers = MAX_PARTLAYERS
  pars_struct%max_directions = MAX_DIRECTIONS
  pars_struct%max_brdf_kernels = MAX_BRDF_KERNELS
  pars_struct%max_brdf_parameters = MAX_BRDF_PARAMETERS
  pars_struct%maxstreams_brdf = MAXSTREAMS_BRDF
  pars_struct%max_msrs_muquad = MAX_MSRS_MUQUAD
  pars_struct%max_msrs_phiquad = MAX_MSRS_PHIQUAD
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
  pars_struct%breonveg_idx = BREONVEG_IDX
  pars_struct%breonsoil_idx = BREONSOIL_IDX
  pars_struct%maxbrdf_idx = MAXBRDF_IDX
  
end subroutine set_lidort_pars



! Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_linsup_inputs_def, only : brdf_linsup_inputs

  type(brdf_linsup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_linsup_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_linsup_inputs_c_alloc_init

! Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
! Initializes only with no allocation
subroutine brdf_linsup_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_linsup_inputs_def
  use lidort_pars

  type(brdf_linsup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  logical(kind=4), dimension(:), pointer :: bs_do_kernel_factor_wfs_ptr
  logical(kind=4), dimension(:,:), pointer :: bs_do_kernel_params_wfs_ptr
  logical(kind=4), dimension(:), pointer :: bs_do_kparams_derivs_ptr
  integer(c_int), pointer :: bs_n_surface_wfs_ptr
  integer(c_int), pointer :: bs_n_kernel_factor_wfs_ptr
  integer(c_int), pointer :: bs_n_kernel_params_wfs_ptr
  

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
  
  
  
end subroutine brdf_linsup_inputs_c_init_only

subroutine brdf_linsup_inputs_c_destroy(fortran_type_c) bind(C)
  use brdf_linsup_inputs_def, only : brdf_linsup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_linsup_inputs_c_destroy

subroutine brdf_linsup_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_linsup_inputs_def, only : brdf_linsup_inputs

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
  

end subroutine brdf_linsup_inputs_c_copy

! Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_linsup_outputs_def, only : brdf_linsup_outputs

  type(brdf_linsup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_linsup_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_linsup_outputs_c_alloc_init

! Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_linsup_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_linsup_outputs_def
  use lidort_pars

  type(brdf_linsup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_brdf_f_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: bs_ls_user_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_ls_emissivity_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_ls_user_emissivity_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_ls_exactdb_brdfunc = 0_fpk
  bs_ls_exactdb_brdfunc_ptr => fortran_type_f%bs_ls_exactdb_brdfunc
  transfer_struct_c%bs_ls_exactdb_brdfunc = c_loc(bs_ls_exactdb_brdfunc_ptr(&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,1),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,2),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,3),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,4)))
  inquire(iolength=transfer_struct_c%bs_ls_exactdb_brdfunc_f_byte_size) fortran_type_f%bs_ls_exactdb_brdfunc(&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,1),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,2),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,3),&
    lbound(fortran_type_f%bs_ls_exactdb_brdfunc,4))
#ifdef ifort
  transfer_struct_c%bs_ls_exactdb_brdfunc_f_byte_size = transfer_struct_c%bs_ls_exactdb_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_exactdb_brdfunc_f_shapes(1) = size(fortran_type_f%bs_ls_exactdb_brdfunc, 1)
  transfer_struct_c%bs_ls_exactdb_brdfunc_f_shapes(2) = size(fortran_type_f%bs_ls_exactdb_brdfunc, 2)
  transfer_struct_c%bs_ls_exactdb_brdfunc_f_shapes(3) = size(fortran_type_f%bs_ls_exactdb_brdfunc, 3)
  transfer_struct_c%bs_ls_exactdb_brdfunc_f_shapes(4) = size(fortran_type_f%bs_ls_exactdb_brdfunc, 4)
  
  
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
    lbound(fortran_type_f%bs_ls_emissivity,2),&
    lbound(fortran_type_f%bs_ls_emissivity,3)))
  inquire(iolength=transfer_struct_c%bs_ls_emissivity_f_byte_size) fortran_type_f%bs_ls_emissivity(&
    lbound(fortran_type_f%bs_ls_emissivity,1),&
    lbound(fortran_type_f%bs_ls_emissivity,2),&
    lbound(fortran_type_f%bs_ls_emissivity,3))
#ifdef ifort
  transfer_struct_c%bs_ls_emissivity_f_byte_size = transfer_struct_c%bs_ls_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_emissivity_f_shapes(1) = size(fortran_type_f%bs_ls_emissivity, 1)
  transfer_struct_c%bs_ls_emissivity_f_shapes(2) = size(fortran_type_f%bs_ls_emissivity, 2)
  transfer_struct_c%bs_ls_emissivity_f_shapes(3) = size(fortran_type_f%bs_ls_emissivity, 3)
  
  
  fortran_type_f%bs_ls_user_emissivity = 0_fpk
  bs_ls_user_emissivity_ptr => fortran_type_f%bs_ls_user_emissivity
  transfer_struct_c%bs_ls_user_emissivity = c_loc(bs_ls_user_emissivity_ptr(&
    lbound(fortran_type_f%bs_ls_user_emissivity,1),&
    lbound(fortran_type_f%bs_ls_user_emissivity,2),&
    lbound(fortran_type_f%bs_ls_user_emissivity,3)))
  inquire(iolength=transfer_struct_c%bs_ls_user_emissivity_f_byte_size) fortran_type_f%bs_ls_user_emissivity(&
    lbound(fortran_type_f%bs_ls_user_emissivity,1),&
    lbound(fortran_type_f%bs_ls_user_emissivity,2),&
    lbound(fortran_type_f%bs_ls_user_emissivity,3))
#ifdef ifort
  transfer_struct_c%bs_ls_user_emissivity_f_byte_size = transfer_struct_c%bs_ls_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_ls_user_emissivity_f_shapes(1) = size(fortran_type_f%bs_ls_user_emissivity, 1)
  transfer_struct_c%bs_ls_user_emissivity_f_shapes(2) = size(fortran_type_f%bs_ls_user_emissivity, 2)
  transfer_struct_c%bs_ls_user_emissivity_f_shapes(3) = size(fortran_type_f%bs_ls_user_emissivity, 3)
  
  
end subroutine brdf_linsup_outputs_c_init_only

subroutine brdf_linsup_outputs_c_destroy(fortran_type_c) bind(C)
  use brdf_linsup_outputs_def, only : brdf_linsup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_linsup_outputs_c_destroy

subroutine brdf_linsup_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_linsup_outputs_def, only : brdf_linsup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_linsup_outputs), pointer :: fortran_type_f_from
  type(brdf_linsup_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_ls_exactdb_brdfunc = fortran_type_f_from%bs_ls_exactdb_brdfunc
  fortran_type_f_to%bs_ls_brdf_f_0 = fortran_type_f_from%bs_ls_brdf_f_0
  fortran_type_f_to%bs_ls_brdf_f = fortran_type_f_from%bs_ls_brdf_f
  fortran_type_f_to%bs_ls_user_brdf_f_0 = fortran_type_f_from%bs_ls_user_brdf_f_0
  fortran_type_f_to%bs_ls_user_brdf_f = fortran_type_f_from%bs_ls_user_brdf_f
  fortran_type_f_to%bs_ls_emissivity = fortran_type_f_from%bs_ls_emissivity
  fortran_type_f_to%bs_ls_user_emissivity = fortran_type_f_from%bs_ls_user_emissivity
  

end subroutine brdf_linsup_outputs_c_copy

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_inputs_def, only : brdf_sup_inputs

  type(brdf_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_sup_inputs_c_alloc_init

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
! Initializes only with no allocation
subroutine brdf_sup_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_inputs_def
  use lidort_pars

  type(brdf_sup_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  logical(kind=4), pointer :: bs_do_user_streams_ptr
  logical(kind=4), pointer :: bs_do_brdf_surface_ptr
  logical(kind=4), pointer :: bs_do_surface_emission_ptr
  integer(c_int), pointer :: bs_nstreams_ptr
  integer(c_int), pointer :: bs_nbeams_ptr
  real(c_double), dimension(:), pointer :: bs_beam_szas_ptr
  integer(c_int), pointer :: bs_n_user_relazms_ptr
  real(c_double), dimension(:), pointer :: bs_user_relazms_ptr
  integer(c_int), pointer :: bs_n_user_streams_ptr
  real(c_double), dimension(:), pointer :: bs_user_angles_input_ptr
  integer(c_int), pointer :: bs_n_brdf_kernels_ptr
  integer(c_int), dimension(:), pointer :: bs_which_brdf_ptr
  integer(c_int), dimension(:), pointer :: bs_n_brdf_parameters_ptr
  real(c_double), dimension(:,:), pointer :: bs_brdf_parameters_ptr
  logical(kind=4), dimension(:), pointer :: bs_lambertian_kernel_flag_ptr
  real(c_double), dimension(:), pointer :: bs_brdf_factors_ptr
  integer(c_int), pointer :: bs_nstreams_brdf_ptr
  logical(kind=4), pointer :: bs_do_shadow_effect_ptr
  logical(kind=4), pointer :: bs_do_exactonly_ptr
  logical(kind=4), pointer :: bs_do_glitter_msrcorr_ptr
  logical(kind=4), pointer :: bs_do_glitter_msrcorr_exactonly_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_order_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_nmuquad_ptr
  integer(c_int), pointer :: bs_glitter_msrcorr_nphiquad_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_do_user_streams = .FALSE.
  bs_do_user_streams_ptr => fortran_type_f%bs_do_user_streams
  transfer_struct_c%bs_do_user_streams = c_loc(bs_do_user_streams_ptr)
  inquire(iolength=transfer_struct_c%bs_do_user_streams_f_byte_size) fortran_type_f%bs_do_user_streams
#ifdef ifort
  transfer_struct_c%bs_do_user_streams_f_byte_size = transfer_struct_c%bs_do_user_streams_f_byte_size * 4
#endif
  
  
  
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
  
  
  
  fortran_type_f%bs_do_exactonly = .FALSE.
  bs_do_exactonly_ptr => fortran_type_f%bs_do_exactonly
  transfer_struct_c%bs_do_exactonly = c_loc(bs_do_exactonly_ptr)
  inquire(iolength=transfer_struct_c%bs_do_exactonly_f_byte_size) fortran_type_f%bs_do_exactonly
#ifdef ifort
  transfer_struct_c%bs_do_exactonly_f_byte_size = transfer_struct_c%bs_do_exactonly_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_glitter_msrcorr = .FALSE.
  bs_do_glitter_msrcorr_ptr => fortran_type_f%bs_do_glitter_msrcorr
  transfer_struct_c%bs_do_glitter_msrcorr = c_loc(bs_do_glitter_msrcorr_ptr)
  inquire(iolength=transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size) fortran_type_f%bs_do_glitter_msrcorr
#ifdef ifort
  transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size = transfer_struct_c%bs_do_glitter_msrcorr_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%bs_do_glitter_msrcorr_exactonly = .FALSE.
  bs_do_glitter_msrcorr_exactonly_ptr => fortran_type_f%bs_do_glitter_msrcorr_exactonly
  transfer_struct_c%bs_do_glitter_msrcorr_exactonly = c_loc(bs_do_glitter_msrcorr_exactonly_ptr)
  inquire(iolength=transfer_struct_c%bs_do_glitter_msrcorr_exactonly_f_byte_size) fortran_type_f%bs_do_glitter_msrcorr_exactonly
#ifdef ifort
  transfer_struct_c%bs_do_glitter_msrcorr_exactonly_f_byte_size = transfer_struct_c%bs_do_glitter_msrcorr_exactonly_f_byte_size * 4
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
  use brdf_sup_inputs_def, only : brdf_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_sup_inputs_c_destroy

subroutine brdf_sup_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_inputs_def, only : brdf_sup_inputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_sup_inputs), pointer :: fortran_type_f_from
  type(brdf_sup_inputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_do_user_streams = fortran_type_f_from%bs_do_user_streams
  fortran_type_f_to%bs_do_brdf_surface = fortran_type_f_from%bs_do_brdf_surface
  fortran_type_f_to%bs_do_surface_emission = fortran_type_f_from%bs_do_surface_emission
  fortran_type_f_to%bs_nstreams = fortran_type_f_from%bs_nstreams
  fortran_type_f_to%bs_nbeams = fortran_type_f_from%bs_nbeams
  fortran_type_f_to%bs_beam_szas = fortran_type_f_from%bs_beam_szas
  fortran_type_f_to%bs_n_user_relazms = fortran_type_f_from%bs_n_user_relazms
  fortran_type_f_to%bs_user_relazms = fortran_type_f_from%bs_user_relazms
  fortran_type_f_to%bs_n_user_streams = fortran_type_f_from%bs_n_user_streams
  fortran_type_f_to%bs_user_angles_input = fortran_type_f_from%bs_user_angles_input
  fortran_type_f_to%bs_n_brdf_kernels = fortran_type_f_from%bs_n_brdf_kernels
  fortran_type_f_to%bs_brdf_names = fortran_type_f_from%bs_brdf_names
  fortran_type_f_to%bs_which_brdf = fortran_type_f_from%bs_which_brdf
  fortran_type_f_to%bs_n_brdf_parameters = fortran_type_f_from%bs_n_brdf_parameters
  fortran_type_f_to%bs_brdf_parameters = fortran_type_f_from%bs_brdf_parameters
  fortran_type_f_to%bs_lambertian_kernel_flag = fortran_type_f_from%bs_lambertian_kernel_flag
  fortran_type_f_to%bs_brdf_factors = fortran_type_f_from%bs_brdf_factors
  fortran_type_f_to%bs_nstreams_brdf = fortran_type_f_from%bs_nstreams_brdf
  fortran_type_f_to%bs_do_shadow_effect = fortran_type_f_from%bs_do_shadow_effect
  fortran_type_f_to%bs_do_exactonly = fortran_type_f_from%bs_do_exactonly
  fortran_type_f_to%bs_do_glitter_msrcorr = fortran_type_f_from%bs_do_glitter_msrcorr
  fortran_type_f_to%bs_do_glitter_msrcorr_exactonly = fortran_type_f_from%bs_do_glitter_msrcorr_exactonly
  fortran_type_f_to%bs_glitter_msrcorr_order = fortran_type_f_from%bs_glitter_msrcorr_order
  fortran_type_f_to%bs_glitter_msrcorr_nmuquad = fortran_type_f_from%bs_glitter_msrcorr_nmuquad
  fortran_type_f_to%bs_glitter_msrcorr_nphiquad = fortran_type_f_from%bs_glitter_msrcorr_nphiquad
  

end subroutine brdf_sup_inputs_c_copy

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_sup_outputs

  type(brdf_sup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_sup_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_sup_outputs_c_alloc_init

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_sup_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def
  use lidort_pars

  type(brdf_sup_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: bs_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: bs_user_brdf_f_ptr
  real(c_double), dimension(:,:), pointer :: bs_emissivity_ptr
  real(c_double), dimension(:,:), pointer :: bs_user_emissivity_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%bs_exactdb_brdfunc = 0_fpk
  bs_exactdb_brdfunc_ptr => fortran_type_f%bs_exactdb_brdfunc
  transfer_struct_c%bs_exactdb_brdfunc = c_loc(bs_exactdb_brdfunc_ptr(&
    lbound(fortran_type_f%bs_exactdb_brdfunc,1),&
    lbound(fortran_type_f%bs_exactdb_brdfunc,2),&
    lbound(fortran_type_f%bs_exactdb_brdfunc,3)))
  inquire(iolength=transfer_struct_c%bs_exactdb_brdfunc_f_byte_size) fortran_type_f%bs_exactdb_brdfunc(&
    lbound(fortran_type_f%bs_exactdb_brdfunc,1),&
    lbound(fortran_type_f%bs_exactdb_brdfunc,2),&
    lbound(fortran_type_f%bs_exactdb_brdfunc,3))
#ifdef ifort
  transfer_struct_c%bs_exactdb_brdfunc_f_byte_size = transfer_struct_c%bs_exactdb_brdfunc_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_exactdb_brdfunc_f_shapes(1) = size(fortran_type_f%bs_exactdb_brdfunc, 1)
  transfer_struct_c%bs_exactdb_brdfunc_f_shapes(2) = size(fortran_type_f%bs_exactdb_brdfunc, 2)
  transfer_struct_c%bs_exactdb_brdfunc_f_shapes(3) = size(fortran_type_f%bs_exactdb_brdfunc, 3)
  
  
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
    lbound(fortran_type_f%bs_emissivity,1),&
    lbound(fortran_type_f%bs_emissivity,2)))
  inquire(iolength=transfer_struct_c%bs_emissivity_f_byte_size) fortran_type_f%bs_emissivity(&
    lbound(fortran_type_f%bs_emissivity,1),&
    lbound(fortran_type_f%bs_emissivity,2))
#ifdef ifort
  transfer_struct_c%bs_emissivity_f_byte_size = transfer_struct_c%bs_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_emissivity_f_shapes(1) = size(fortran_type_f%bs_emissivity, 1)
  transfer_struct_c%bs_emissivity_f_shapes(2) = size(fortran_type_f%bs_emissivity, 2)
  
  
  fortran_type_f%bs_user_emissivity = 0_fpk
  bs_user_emissivity_ptr => fortran_type_f%bs_user_emissivity
  transfer_struct_c%bs_user_emissivity = c_loc(bs_user_emissivity_ptr(&
    lbound(fortran_type_f%bs_user_emissivity,1),&
    lbound(fortran_type_f%bs_user_emissivity,2)))
  inquire(iolength=transfer_struct_c%bs_user_emissivity_f_byte_size) fortran_type_f%bs_user_emissivity(&
    lbound(fortran_type_f%bs_user_emissivity,1),&
    lbound(fortran_type_f%bs_user_emissivity,2))
#ifdef ifort
  transfer_struct_c%bs_user_emissivity_f_byte_size = transfer_struct_c%bs_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%bs_user_emissivity_f_shapes(1) = size(fortran_type_f%bs_user_emissivity, 1)
  transfer_struct_c%bs_user_emissivity_f_shapes(2) = size(fortran_type_f%bs_user_emissivity, 2)
  
  
end subroutine brdf_sup_outputs_c_init_only

subroutine brdf_sup_outputs_c_destroy(fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_sup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_sup_outputs_c_destroy

subroutine brdf_sup_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_outputs_def, only : brdf_sup_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(brdf_sup_outputs), pointer :: fortran_type_f_from
  type(brdf_sup_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%bs_exactdb_brdfunc = fortran_type_f_from%bs_exactdb_brdfunc
  fortran_type_f_to%bs_brdf_f_0 = fortran_type_f_from%bs_brdf_f_0
  fortran_type_f_to%bs_brdf_f = fortran_type_f_from%bs_brdf_f
  fortran_type_f_to%bs_user_brdf_f_0 = fortran_type_f_from%bs_user_brdf_f_0
  fortran_type_f_to%bs_user_brdf_f = fortran_type_f_from%bs_user_brdf_f
  fortran_type_f_to%bs_emissivity = fortran_type_f_from%bs_emissivity
  fortran_type_f_to%bs_user_emissivity = fortran_type_f_from%bs_user_emissivity
  

end subroutine brdf_sup_outputs_c_copy

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_input_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

  type(brdf_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call brdf_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine brdf_input_exception_handling_c_alloc_init

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Initializes only with no allocation
subroutine brdf_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use brdf_sup_outputs_def
  use lidort_pars

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
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine brdf_input_exception_handling_c_destroy

subroutine brdf_input_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

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

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lincontrol_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lincontrol

  type(lidort_fixed_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_lincontrol_c_alloc_init

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_lincontrol_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def
  use lidort_pars

  type(lidort_fixed_lincontrol_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_column_linearization_ptr
  logical(kind=4), pointer :: ts_do_profile_linearization_ptr
  logical(kind=4), pointer :: ts_do_surface_linearization_ptr
  logical(kind=4), pointer :: ts_do_sleave_wfs_ptr
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
  
  
  
  fortran_type_f%ts_do_surface_linearization = .FALSE.
  ts_do_surface_linearization_ptr => fortran_type_f%ts_do_surface_linearization
  transfer_struct_c%ts_do_surface_linearization = c_loc(ts_do_surface_linearization_ptr)
  inquire(iolength=transfer_struct_c%ts_do_surface_linearization_f_byte_size) fortran_type_f%ts_do_surface_linearization
#ifdef ifort
  transfer_struct_c%ts_do_surface_linearization_f_byte_size = transfer_struct_c%ts_do_surface_linearization_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sleave_wfs = .FALSE.
  ts_do_sleave_wfs_ptr => fortran_type_f%ts_do_sleave_wfs
  transfer_struct_c%ts_do_sleave_wfs = c_loc(ts_do_sleave_wfs_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sleave_wfs_f_byte_size) fortran_type_f%ts_do_sleave_wfs
#ifdef ifort
  transfer_struct_c%ts_do_sleave_wfs_f_byte_size = transfer_struct_c%ts_do_sleave_wfs_f_byte_size * 4
#endif
  
  
  
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
  
  
  
end subroutine lidort_fixed_lincontrol_c_init_only

subroutine lidort_fixed_lincontrol_c_destroy(fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_lincontrol_c_destroy

subroutine lidort_fixed_lincontrol_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lincontrol

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f_from
  type(lidort_fixed_lincontrol), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_column_linearization = fortran_type_f_from%ts_do_column_linearization
  fortran_type_f_to%ts_do_profile_linearization = fortran_type_f_from%ts_do_profile_linearization
  fortran_type_f_to%ts_do_surface_linearization = fortran_type_f_from%ts_do_surface_linearization
  fortran_type_f_to%ts_do_sleave_wfs = fortran_type_f_from%ts_do_sleave_wfs
  fortran_type_f_to%ts_layer_vary_flag = fortran_type_f_from%ts_layer_vary_flag
  fortran_type_f_to%ts_layer_vary_number = fortran_type_f_from%ts_layer_vary_number
  fortran_type_f_to%ts_n_totalcolumn_wfs = fortran_type_f_from%ts_n_totalcolumn_wfs
  fortran_type_f_to%ts_n_surface_wfs = fortran_type_f_from%ts_n_surface_wfs
  fortran_type_f_to%ts_n_sleave_wfs = fortran_type_f_from%ts_n_sleave_wfs
  

end subroutine lidort_fixed_lincontrol_c_copy

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_linoptical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_linoptical

  type(lidort_fixed_linoptical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_linoptical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_linoptical_c_alloc_init

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_linoptical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def
  use lidort_pars

  type(lidort_fixed_linoptical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_l_deltau_vert_input_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_l_omega_total_input_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_l_phasmoms_total_input_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_l_deltau_vert_input = 0_fpk
  ts_l_deltau_vert_input_ptr => fortran_type_f%ts_l_deltau_vert_input
  transfer_struct_c%ts_l_deltau_vert_input = c_loc(ts_l_deltau_vert_input_ptr(&
    lbound(fortran_type_f%ts_l_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,2),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,3)))
  inquire(iolength=transfer_struct_c%ts_l_deltau_vert_input_f_byte_size) fortran_type_f%ts_l_deltau_vert_input(&
    lbound(fortran_type_f%ts_l_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,2),&
    lbound(fortran_type_f%ts_l_deltau_vert_input,3))
#ifdef ifort
  transfer_struct_c%ts_l_deltau_vert_input_f_byte_size = transfer_struct_c%ts_l_deltau_vert_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_deltau_vert_input_f_shapes(1) = size(fortran_type_f%ts_l_deltau_vert_input, 1)
  transfer_struct_c%ts_l_deltau_vert_input_f_shapes(2) = size(fortran_type_f%ts_l_deltau_vert_input, 2)
  transfer_struct_c%ts_l_deltau_vert_input_f_shapes(3) = size(fortran_type_f%ts_l_deltau_vert_input, 3)
  
  
  fortran_type_f%ts_l_omega_total_input = 0_fpk
  ts_l_omega_total_input_ptr => fortran_type_f%ts_l_omega_total_input
  transfer_struct_c%ts_l_omega_total_input = c_loc(ts_l_omega_total_input_ptr(&
    lbound(fortran_type_f%ts_l_omega_total_input,1),&
    lbound(fortran_type_f%ts_l_omega_total_input,2),&
    lbound(fortran_type_f%ts_l_omega_total_input,3)))
  inquire(iolength=transfer_struct_c%ts_l_omega_total_input_f_byte_size) fortran_type_f%ts_l_omega_total_input(&
    lbound(fortran_type_f%ts_l_omega_total_input,1),&
    lbound(fortran_type_f%ts_l_omega_total_input,2),&
    lbound(fortran_type_f%ts_l_omega_total_input,3))
#ifdef ifort
  transfer_struct_c%ts_l_omega_total_input_f_byte_size = transfer_struct_c%ts_l_omega_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_omega_total_input_f_shapes(1) = size(fortran_type_f%ts_l_omega_total_input, 1)
  transfer_struct_c%ts_l_omega_total_input_f_shapes(2) = size(fortran_type_f%ts_l_omega_total_input, 2)
  transfer_struct_c%ts_l_omega_total_input_f_shapes(3) = size(fortran_type_f%ts_l_omega_total_input, 3)
  
  
  fortran_type_f%ts_l_phasmoms_total_input = 0_fpk
  ts_l_phasmoms_total_input_ptr => fortran_type_f%ts_l_phasmoms_total_input
  transfer_struct_c%ts_l_phasmoms_total_input = c_loc(ts_l_phasmoms_total_input_ptr(&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,3),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,4)))
  inquire(iolength=transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size) fortran_type_f%ts_l_phasmoms_total_input(&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,3),&
    lbound(fortran_type_f%ts_l_phasmoms_total_input,4))
#ifdef ifort
  transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size = transfer_struct_c%ts_l_phasmoms_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(1) = size(fortran_type_f%ts_l_phasmoms_total_input, 1)
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(2) = size(fortran_type_f%ts_l_phasmoms_total_input, 2)
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(3) = size(fortran_type_f%ts_l_phasmoms_total_input, 3)
  transfer_struct_c%ts_l_phasmoms_total_input_f_shapes(4) = size(fortran_type_f%ts_l_phasmoms_total_input, 4)
  
  
end subroutine lidort_fixed_linoptical_c_init_only

subroutine lidort_fixed_linoptical_c_destroy(fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_linoptical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_linoptical_c_destroy

subroutine lidort_fixed_linoptical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_linoptical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_linoptical), pointer :: fortran_type_f_from
  type(lidort_fixed_linoptical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_l_deltau_vert_input = fortran_type_f_from%ts_l_deltau_vert_input
  fortran_type_f_to%ts_l_omega_total_input = fortran_type_f_from%ts_l_omega_total_input
  fortran_type_f_to%ts_l_phasmoms_total_input = fortran_type_f_from%ts_l_phasmoms_total_input
  

end subroutine lidort_fixed_linoptical_c_copy

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lininputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lininputs

  type(lidort_fixed_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_lininputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_lininputs_c_alloc_init

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_lininputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def
  use lidort_pars

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
  use lidort_lininputs_def, only : lidort_fixed_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_lininputs_c_destroy

subroutine lidort_fixed_lininputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_lininputs), pointer :: fortran_type_f_from
  type(lidort_fixed_lininputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%cont = fortran_type_f_from%cont
  fortran_type_f_to%optical = fortran_type_f_from%optical
  

end subroutine lidort_fixed_lininputs_c_copy

! Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_lininputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_modified_lininputs

  type(lidort_modified_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_lininputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_lininputs_c_alloc_init

! Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_lininputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_lininputs_def
  use lidort_pars

  type(lidort_modified_lininputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  integer(c_int), pointer :: dummy_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%dummy = 0
  dummy_ptr => fortran_type_f%dummy
  transfer_struct_c%dummy = c_loc(dummy_ptr)
  inquire(iolength=transfer_struct_c%dummy_f_byte_size) fortran_type_f%dummy
#ifdef ifort
  transfer_struct_c%dummy_f_byte_size = transfer_struct_c%dummy_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_modified_lininputs_c_init_only

subroutine lidort_modified_lininputs_c_destroy(fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_modified_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_lininputs_c_destroy

subroutine lidort_modified_lininputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_lininputs_def, only : lidort_modified_lininputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_lininputs), pointer :: fortran_type_f_from
  type(lidort_modified_lininputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%dummy = fortran_type_f_from%dummy
  

end subroutine lidort_modified_lininputs_c_copy

! Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linatmos_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linatmos

  type(lidort_linatmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linatmos_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linatmos_c_alloc_init

! Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linatmos_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def
  use lidort_pars

  type(lidort_linatmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:,:), pointer :: ts_columnwf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_mint_columnwf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_flux_columnwf_ptr
  real(c_double), dimension(:,:,:,:,:,:), pointer :: ts_profilewf_ptr
  real(c_double), dimension(:,:,:,:,:,:), pointer :: ts_mint_profilewf_ptr
  real(c_double), dimension(:,:,:,:,:,:), pointer :: ts_flux_profilewf_ptr
  

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
    lbound(fortran_type_f%ts_columnwf,4),&
    lbound(fortran_type_f%ts_columnwf,5)))
  inquire(iolength=transfer_struct_c%ts_columnwf_f_byte_size) fortran_type_f%ts_columnwf(&
    lbound(fortran_type_f%ts_columnwf,1),&
    lbound(fortran_type_f%ts_columnwf,2),&
    lbound(fortran_type_f%ts_columnwf,3),&
    lbound(fortran_type_f%ts_columnwf,4),&
    lbound(fortran_type_f%ts_columnwf,5))
#ifdef ifort
  transfer_struct_c%ts_columnwf_f_byte_size = transfer_struct_c%ts_columnwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_columnwf_f_shapes(1) = size(fortran_type_f%ts_columnwf, 1)
  transfer_struct_c%ts_columnwf_f_shapes(2) = size(fortran_type_f%ts_columnwf, 2)
  transfer_struct_c%ts_columnwf_f_shapes(3) = size(fortran_type_f%ts_columnwf, 3)
  transfer_struct_c%ts_columnwf_f_shapes(4) = size(fortran_type_f%ts_columnwf, 4)
  transfer_struct_c%ts_columnwf_f_shapes(5) = size(fortran_type_f%ts_columnwf, 5)
  
  
  fortran_type_f%ts_mint_columnwf = 0_fpk
  ts_mint_columnwf_ptr => fortran_type_f%ts_mint_columnwf
  transfer_struct_c%ts_mint_columnwf = c_loc(ts_mint_columnwf_ptr(&
    lbound(fortran_type_f%ts_mint_columnwf,1),&
    lbound(fortran_type_f%ts_mint_columnwf,2),&
    lbound(fortran_type_f%ts_mint_columnwf,3),&
    lbound(fortran_type_f%ts_mint_columnwf,4),&
    lbound(fortran_type_f%ts_mint_columnwf,5)))
  inquire(iolength=transfer_struct_c%ts_mint_columnwf_f_byte_size) fortran_type_f%ts_mint_columnwf(&
    lbound(fortran_type_f%ts_mint_columnwf,1),&
    lbound(fortran_type_f%ts_mint_columnwf,2),&
    lbound(fortran_type_f%ts_mint_columnwf,3),&
    lbound(fortran_type_f%ts_mint_columnwf,4),&
    lbound(fortran_type_f%ts_mint_columnwf,5))
#ifdef ifort
  transfer_struct_c%ts_mint_columnwf_f_byte_size = transfer_struct_c%ts_mint_columnwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_mint_columnwf_f_shapes(1) = size(fortran_type_f%ts_mint_columnwf, 1)
  transfer_struct_c%ts_mint_columnwf_f_shapes(2) = size(fortran_type_f%ts_mint_columnwf, 2)
  transfer_struct_c%ts_mint_columnwf_f_shapes(3) = size(fortran_type_f%ts_mint_columnwf, 3)
  transfer_struct_c%ts_mint_columnwf_f_shapes(4) = size(fortran_type_f%ts_mint_columnwf, 4)
  transfer_struct_c%ts_mint_columnwf_f_shapes(5) = size(fortran_type_f%ts_mint_columnwf, 5)
  
  
  fortran_type_f%ts_flux_columnwf = 0_fpk
  ts_flux_columnwf_ptr => fortran_type_f%ts_flux_columnwf
  transfer_struct_c%ts_flux_columnwf = c_loc(ts_flux_columnwf_ptr(&
    lbound(fortran_type_f%ts_flux_columnwf,1),&
    lbound(fortran_type_f%ts_flux_columnwf,2),&
    lbound(fortran_type_f%ts_flux_columnwf,3),&
    lbound(fortran_type_f%ts_flux_columnwf,4),&
    lbound(fortran_type_f%ts_flux_columnwf,5)))
  inquire(iolength=transfer_struct_c%ts_flux_columnwf_f_byte_size) fortran_type_f%ts_flux_columnwf(&
    lbound(fortran_type_f%ts_flux_columnwf,1),&
    lbound(fortran_type_f%ts_flux_columnwf,2),&
    lbound(fortran_type_f%ts_flux_columnwf,3),&
    lbound(fortran_type_f%ts_flux_columnwf,4),&
    lbound(fortran_type_f%ts_flux_columnwf,5))
#ifdef ifort
  transfer_struct_c%ts_flux_columnwf_f_byte_size = transfer_struct_c%ts_flux_columnwf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_columnwf_f_shapes(1) = size(fortran_type_f%ts_flux_columnwf, 1)
  transfer_struct_c%ts_flux_columnwf_f_shapes(2) = size(fortran_type_f%ts_flux_columnwf, 2)
  transfer_struct_c%ts_flux_columnwf_f_shapes(3) = size(fortran_type_f%ts_flux_columnwf, 3)
  transfer_struct_c%ts_flux_columnwf_f_shapes(4) = size(fortran_type_f%ts_flux_columnwf, 4)
  transfer_struct_c%ts_flux_columnwf_f_shapes(5) = size(fortran_type_f%ts_flux_columnwf, 5)
  
  
  fortran_type_f%ts_profilewf = 0_fpk
  ts_profilewf_ptr => fortran_type_f%ts_profilewf
  transfer_struct_c%ts_profilewf = c_loc(ts_profilewf_ptr(&
    lbound(fortran_type_f%ts_profilewf,1),&
    lbound(fortran_type_f%ts_profilewf,2),&
    lbound(fortran_type_f%ts_profilewf,3),&
    lbound(fortran_type_f%ts_profilewf,4),&
    lbound(fortran_type_f%ts_profilewf,5),&
    lbound(fortran_type_f%ts_profilewf,6)))
  inquire(iolength=transfer_struct_c%ts_profilewf_f_byte_size) fortran_type_f%ts_profilewf(&
    lbound(fortran_type_f%ts_profilewf,1),&
    lbound(fortran_type_f%ts_profilewf,2),&
    lbound(fortran_type_f%ts_profilewf,3),&
    lbound(fortran_type_f%ts_profilewf,4),&
    lbound(fortran_type_f%ts_profilewf,5),&
    lbound(fortran_type_f%ts_profilewf,6))
#ifdef ifort
  transfer_struct_c%ts_profilewf_f_byte_size = transfer_struct_c%ts_profilewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_profilewf_f_shapes(1) = size(fortran_type_f%ts_profilewf, 1)
  transfer_struct_c%ts_profilewf_f_shapes(2) = size(fortran_type_f%ts_profilewf, 2)
  transfer_struct_c%ts_profilewf_f_shapes(3) = size(fortran_type_f%ts_profilewf, 3)
  transfer_struct_c%ts_profilewf_f_shapes(4) = size(fortran_type_f%ts_profilewf, 4)
  transfer_struct_c%ts_profilewf_f_shapes(5) = size(fortran_type_f%ts_profilewf, 5)
  transfer_struct_c%ts_profilewf_f_shapes(6) = size(fortran_type_f%ts_profilewf, 6)
  
  
  fortran_type_f%ts_mint_profilewf = 0_fpk
  ts_mint_profilewf_ptr => fortran_type_f%ts_mint_profilewf
  transfer_struct_c%ts_mint_profilewf = c_loc(ts_mint_profilewf_ptr(&
    lbound(fortran_type_f%ts_mint_profilewf,1),&
    lbound(fortran_type_f%ts_mint_profilewf,2),&
    lbound(fortran_type_f%ts_mint_profilewf,3),&
    lbound(fortran_type_f%ts_mint_profilewf,4),&
    lbound(fortran_type_f%ts_mint_profilewf,5),&
    lbound(fortran_type_f%ts_mint_profilewf,6)))
  inquire(iolength=transfer_struct_c%ts_mint_profilewf_f_byte_size) fortran_type_f%ts_mint_profilewf(&
    lbound(fortran_type_f%ts_mint_profilewf,1),&
    lbound(fortran_type_f%ts_mint_profilewf,2),&
    lbound(fortran_type_f%ts_mint_profilewf,3),&
    lbound(fortran_type_f%ts_mint_profilewf,4),&
    lbound(fortran_type_f%ts_mint_profilewf,5),&
    lbound(fortran_type_f%ts_mint_profilewf,6))
#ifdef ifort
  transfer_struct_c%ts_mint_profilewf_f_byte_size = transfer_struct_c%ts_mint_profilewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_mint_profilewf_f_shapes(1) = size(fortran_type_f%ts_mint_profilewf, 1)
  transfer_struct_c%ts_mint_profilewf_f_shapes(2) = size(fortran_type_f%ts_mint_profilewf, 2)
  transfer_struct_c%ts_mint_profilewf_f_shapes(3) = size(fortran_type_f%ts_mint_profilewf, 3)
  transfer_struct_c%ts_mint_profilewf_f_shapes(4) = size(fortran_type_f%ts_mint_profilewf, 4)
  transfer_struct_c%ts_mint_profilewf_f_shapes(5) = size(fortran_type_f%ts_mint_profilewf, 5)
  transfer_struct_c%ts_mint_profilewf_f_shapes(6) = size(fortran_type_f%ts_mint_profilewf, 6)
  
  
  fortran_type_f%ts_flux_profilewf = 0_fpk
  ts_flux_profilewf_ptr => fortran_type_f%ts_flux_profilewf
  transfer_struct_c%ts_flux_profilewf = c_loc(ts_flux_profilewf_ptr(&
    lbound(fortran_type_f%ts_flux_profilewf,1),&
    lbound(fortran_type_f%ts_flux_profilewf,2),&
    lbound(fortran_type_f%ts_flux_profilewf,3),&
    lbound(fortran_type_f%ts_flux_profilewf,4),&
    lbound(fortran_type_f%ts_flux_profilewf,5),&
    lbound(fortran_type_f%ts_flux_profilewf,6)))
  inquire(iolength=transfer_struct_c%ts_flux_profilewf_f_byte_size) fortran_type_f%ts_flux_profilewf(&
    lbound(fortran_type_f%ts_flux_profilewf,1),&
    lbound(fortran_type_f%ts_flux_profilewf,2),&
    lbound(fortran_type_f%ts_flux_profilewf,3),&
    lbound(fortran_type_f%ts_flux_profilewf,4),&
    lbound(fortran_type_f%ts_flux_profilewf,5),&
    lbound(fortran_type_f%ts_flux_profilewf,6))
#ifdef ifort
  transfer_struct_c%ts_flux_profilewf_f_byte_size = transfer_struct_c%ts_flux_profilewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_profilewf_f_shapes(1) = size(fortran_type_f%ts_flux_profilewf, 1)
  transfer_struct_c%ts_flux_profilewf_f_shapes(2) = size(fortran_type_f%ts_flux_profilewf, 2)
  transfer_struct_c%ts_flux_profilewf_f_shapes(3) = size(fortran_type_f%ts_flux_profilewf, 3)
  transfer_struct_c%ts_flux_profilewf_f_shapes(4) = size(fortran_type_f%ts_flux_profilewf, 4)
  transfer_struct_c%ts_flux_profilewf_f_shapes(5) = size(fortran_type_f%ts_flux_profilewf, 5)
  transfer_struct_c%ts_flux_profilewf_f_shapes(6) = size(fortran_type_f%ts_flux_profilewf, 6)
  
  
end subroutine lidort_linatmos_c_init_only

subroutine lidort_linatmos_c_destroy(fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linatmos

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linatmos_c_destroy

subroutine lidort_linatmos_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linoutputs_def, only : lidort_linatmos

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linatmos), pointer :: fortran_type_f_from
  type(lidort_linatmos), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_columnwf = fortran_type_f_from%ts_columnwf
  fortran_type_f_to%ts_mint_columnwf = fortran_type_f_from%ts_mint_columnwf
  fortran_type_f_to%ts_flux_columnwf = fortran_type_f_from%ts_flux_columnwf
  fortran_type_f_to%ts_profilewf = fortran_type_f_from%ts_profilewf
  fortran_type_f_to%ts_mint_profilewf = fortran_type_f_from%ts_mint_profilewf
  fortran_type_f_to%ts_flux_profilewf = fortran_type_f_from%ts_flux_profilewf
  

end subroutine lidort_linatmos_c_copy

! Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linsurf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linsurf

  type(lidort_linsurf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsurf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsurf_c_alloc_init

! Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linsurf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def
  use lidort_pars

  type(lidort_linsurf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:,:), pointer :: ts_surfacewf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_mint_surfacewf_ptr
  real(c_double), dimension(:,:,:,:,:), pointer :: ts_flux_surfacewf_ptr
  

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
    lbound(fortran_type_f%ts_surfacewf,4),&
    lbound(fortran_type_f%ts_surfacewf,5)))
  inquire(iolength=transfer_struct_c%ts_surfacewf_f_byte_size) fortran_type_f%ts_surfacewf(&
    lbound(fortran_type_f%ts_surfacewf,1),&
    lbound(fortran_type_f%ts_surfacewf,2),&
    lbound(fortran_type_f%ts_surfacewf,3),&
    lbound(fortran_type_f%ts_surfacewf,4),&
    lbound(fortran_type_f%ts_surfacewf,5))
#ifdef ifort
  transfer_struct_c%ts_surfacewf_f_byte_size = transfer_struct_c%ts_surfacewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_surfacewf_f_shapes(1) = size(fortran_type_f%ts_surfacewf, 1)
  transfer_struct_c%ts_surfacewf_f_shapes(2) = size(fortran_type_f%ts_surfacewf, 2)
  transfer_struct_c%ts_surfacewf_f_shapes(3) = size(fortran_type_f%ts_surfacewf, 3)
  transfer_struct_c%ts_surfacewf_f_shapes(4) = size(fortran_type_f%ts_surfacewf, 4)
  transfer_struct_c%ts_surfacewf_f_shapes(5) = size(fortran_type_f%ts_surfacewf, 5)
  
  
  fortran_type_f%ts_mint_surfacewf = 0_fpk
  ts_mint_surfacewf_ptr => fortran_type_f%ts_mint_surfacewf
  transfer_struct_c%ts_mint_surfacewf = c_loc(ts_mint_surfacewf_ptr(&
    lbound(fortran_type_f%ts_mint_surfacewf,1),&
    lbound(fortran_type_f%ts_mint_surfacewf,2),&
    lbound(fortran_type_f%ts_mint_surfacewf,3),&
    lbound(fortran_type_f%ts_mint_surfacewf,4),&
    lbound(fortran_type_f%ts_mint_surfacewf,5)))
  inquire(iolength=transfer_struct_c%ts_mint_surfacewf_f_byte_size) fortran_type_f%ts_mint_surfacewf(&
    lbound(fortran_type_f%ts_mint_surfacewf,1),&
    lbound(fortran_type_f%ts_mint_surfacewf,2),&
    lbound(fortran_type_f%ts_mint_surfacewf,3),&
    lbound(fortran_type_f%ts_mint_surfacewf,4),&
    lbound(fortran_type_f%ts_mint_surfacewf,5))
#ifdef ifort
  transfer_struct_c%ts_mint_surfacewf_f_byte_size = transfer_struct_c%ts_mint_surfacewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_mint_surfacewf_f_shapes(1) = size(fortran_type_f%ts_mint_surfacewf, 1)
  transfer_struct_c%ts_mint_surfacewf_f_shapes(2) = size(fortran_type_f%ts_mint_surfacewf, 2)
  transfer_struct_c%ts_mint_surfacewf_f_shapes(3) = size(fortran_type_f%ts_mint_surfacewf, 3)
  transfer_struct_c%ts_mint_surfacewf_f_shapes(4) = size(fortran_type_f%ts_mint_surfacewf, 4)
  transfer_struct_c%ts_mint_surfacewf_f_shapes(5) = size(fortran_type_f%ts_mint_surfacewf, 5)
  
  
  fortran_type_f%ts_flux_surfacewf = 0_fpk
  ts_flux_surfacewf_ptr => fortran_type_f%ts_flux_surfacewf
  transfer_struct_c%ts_flux_surfacewf = c_loc(ts_flux_surfacewf_ptr(&
    lbound(fortran_type_f%ts_flux_surfacewf,1),&
    lbound(fortran_type_f%ts_flux_surfacewf,2),&
    lbound(fortran_type_f%ts_flux_surfacewf,3),&
    lbound(fortran_type_f%ts_flux_surfacewf,4),&
    lbound(fortran_type_f%ts_flux_surfacewf,5)))
  inquire(iolength=transfer_struct_c%ts_flux_surfacewf_f_byte_size) fortran_type_f%ts_flux_surfacewf(&
    lbound(fortran_type_f%ts_flux_surfacewf,1),&
    lbound(fortran_type_f%ts_flux_surfacewf,2),&
    lbound(fortran_type_f%ts_flux_surfacewf,3),&
    lbound(fortran_type_f%ts_flux_surfacewf,4),&
    lbound(fortran_type_f%ts_flux_surfacewf,5))
#ifdef ifort
  transfer_struct_c%ts_flux_surfacewf_f_byte_size = transfer_struct_c%ts_flux_surfacewf_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_surfacewf_f_shapes(1) = size(fortran_type_f%ts_flux_surfacewf, 1)
  transfer_struct_c%ts_flux_surfacewf_f_shapes(2) = size(fortran_type_f%ts_flux_surfacewf, 2)
  transfer_struct_c%ts_flux_surfacewf_f_shapes(3) = size(fortran_type_f%ts_flux_surfacewf, 3)
  transfer_struct_c%ts_flux_surfacewf_f_shapes(4) = size(fortran_type_f%ts_flux_surfacewf, 4)
  transfer_struct_c%ts_flux_surfacewf_f_shapes(5) = size(fortran_type_f%ts_flux_surfacewf, 5)
  
  
end subroutine lidort_linsurf_c_init_only

subroutine lidort_linsurf_c_destroy(fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linsurf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsurf_c_destroy

subroutine lidort_linsurf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linoutputs_def, only : lidort_linsurf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsurf), pointer :: fortran_type_f_from
  type(lidort_linsurf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_surfacewf = fortran_type_f_from%ts_surfacewf
  fortran_type_f_to%ts_mint_surfacewf = fortran_type_f_from%ts_mint_surfacewf
  fortran_type_f_to%ts_flux_surfacewf = fortran_type_f_from%ts_flux_surfacewf
  

end subroutine lidort_linsurf_c_copy

! Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linoutputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linoutputs

  type(lidort_linoutputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linoutputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linoutputs_c_alloc_init

! Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_linoutputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linoutputs_def
  use lidort_pars

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
  use lidort_linoutputs_def, only : lidort_linoutputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linoutputs_c_destroy

subroutine lidort_linoutputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linoutputs_def, only : lidort_linoutputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linoutputs), pointer :: fortran_type_f_from
  type(lidort_linoutputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%atmos = fortran_type_f_from%atmos
  fortran_type_f_to%surf = fortran_type_f_from%surf
  

end subroutine lidort_linoutputs_c_copy

! Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_brdf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_brdf_def, only : lidort_linsup_brdf

  type(lidort_linsup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_brdf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_brdf_c_alloc_init

! Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_brdf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_brdf_def
  use lidort_pars

  type(lidort_linsup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_brdf_f_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_ls_user_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_ls_emissivity_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_ls_user_emissivity_ptr
  

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
    lbound(fortran_type_f%ts_ls_emissivity,2),&
    lbound(fortran_type_f%ts_ls_emissivity,3)))
  inquire(iolength=transfer_struct_c%ts_ls_emissivity_f_byte_size) fortran_type_f%ts_ls_emissivity(&
    lbound(fortran_type_f%ts_ls_emissivity,1),&
    lbound(fortran_type_f%ts_ls_emissivity,2),&
    lbound(fortran_type_f%ts_ls_emissivity,3))
#ifdef ifort
  transfer_struct_c%ts_ls_emissivity_f_byte_size = transfer_struct_c%ts_ls_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_emissivity_f_shapes(1) = size(fortran_type_f%ts_ls_emissivity, 1)
  transfer_struct_c%ts_ls_emissivity_f_shapes(2) = size(fortran_type_f%ts_ls_emissivity, 2)
  transfer_struct_c%ts_ls_emissivity_f_shapes(3) = size(fortran_type_f%ts_ls_emissivity, 3)
  
  
  fortran_type_f%ts_ls_user_emissivity = 0_fpk
  ts_ls_user_emissivity_ptr => fortran_type_f%ts_ls_user_emissivity
  transfer_struct_c%ts_ls_user_emissivity = c_loc(ts_ls_user_emissivity_ptr(&
    lbound(fortran_type_f%ts_ls_user_emissivity,1),&
    lbound(fortran_type_f%ts_ls_user_emissivity,2),&
    lbound(fortran_type_f%ts_ls_user_emissivity,3)))
  inquire(iolength=transfer_struct_c%ts_ls_user_emissivity_f_byte_size) fortran_type_f%ts_ls_user_emissivity(&
    lbound(fortran_type_f%ts_ls_user_emissivity,1),&
    lbound(fortran_type_f%ts_ls_user_emissivity,2),&
    lbound(fortran_type_f%ts_ls_user_emissivity,3))
#ifdef ifort
  transfer_struct_c%ts_ls_user_emissivity_f_byte_size = transfer_struct_c%ts_ls_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_ls_user_emissivity_f_shapes(1) = size(fortran_type_f%ts_ls_user_emissivity, 1)
  transfer_struct_c%ts_ls_user_emissivity_f_shapes(2) = size(fortran_type_f%ts_ls_user_emissivity, 2)
  transfer_struct_c%ts_ls_user_emissivity_f_shapes(3) = size(fortran_type_f%ts_ls_user_emissivity, 3)
  
  
end subroutine lidort_linsup_brdf_c_init_only

subroutine lidort_linsup_brdf_c_destroy(fortran_type_c) bind(C)
  use lidort_linsup_brdf_def, only : lidort_linsup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_brdf_c_destroy

subroutine lidort_linsup_brdf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_brdf_def, only : lidort_linsup_brdf

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

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_atmos_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_atmos

  type(lidort_linsup_ss_atmos_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_atmos_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_atmos_c_alloc_init

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_atmos_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def
  use lidort_pars

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
  use lidort_linsup_ss_def, only : lidort_linsup_ss_atmos

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_atmos_c_destroy

subroutine lidort_linsup_ss_atmos_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_atmos

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

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_surf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_surf

  type(lidort_linsup_ss_surf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_surf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_surf_c_alloc_init

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_surf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def
  use lidort_pars

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
  use lidort_linsup_ss_def, only : lidort_linsup_ss_surf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_surf_c_destroy

subroutine lidort_linsup_ss_surf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_surf

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f_from
  type(lidort_linsup_ss_surf), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_surfacewf_db = fortran_type_f_from%ts_surfacewf_db
  

end subroutine lidort_linsup_ss_surf_c_copy

! Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss

  type(lidort_linsup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_ss_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_ss_c_alloc_init

! Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_ss_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_ss_def
  use lidort_pars

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
  use lidort_linsup_ss_def, only : lidort_linsup_ss

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_ss_c_destroy

subroutine lidort_linsup_ss_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_linsup_ss), pointer :: fortran_type_f_from
  type(lidort_linsup_ss), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%atmos = fortran_type_f_from%atmos
  fortran_type_f_to%surf = fortran_type_f_from%surf
  

end subroutine lidort_linsup_ss_c_copy

! Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_sleave_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_sleave_def, only : lidort_linsup_sleave

  type(lidort_linsup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_sleave_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_sleave_c_alloc_init

! Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_sleave_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_sleave_def
  use lidort_pars

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
  use lidort_linsup_sleave_def, only : lidort_linsup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_sleave_c_destroy

subroutine lidort_linsup_sleave_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_sleave_def, only : lidort_linsup_sleave

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

! Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_inout_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_inout_def, only : lidort_linsup_inout

  type(lidort_linsup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_linsup_inout_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_linsup_inout_c_alloc_init

! Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
! Initializes only with no allocation
subroutine lidort_linsup_inout_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_linsup_inout_def
  use lidort_linsup_brdf_def
  use lidort_linsup_ss_def
  use lidort_linsup_sleave_def

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
  use lidort_linsup_inout_def, only : lidort_linsup_inout

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_linsup_inout_c_destroy

subroutine lidort_linsup_inout_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_linsup_inout_def, only : lidort_linsup_inout

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

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_main_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_main_outputs

  type(lidort_main_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_main_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_main_outputs_c_alloc_init

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_main_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def
  use lidort_pars

  type(lidort_main_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:,:), pointer :: ts_intensity_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_mean_intensity_ptr
  real(c_double), dimension(:,:,:,:), pointer :: ts_flux_integral_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_dnflux_direct_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_dnmean_direct_ptr
  integer(c_int), dimension(:,:), pointer :: ts_fourier_saved_ptr
  integer(c_int), pointer :: ts_n_geometries_ptr
  

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
    lbound(fortran_type_f%ts_intensity,3),&
    lbound(fortran_type_f%ts_intensity,4)))
  inquire(iolength=transfer_struct_c%ts_intensity_f_byte_size) fortran_type_f%ts_intensity(&
    lbound(fortran_type_f%ts_intensity,1),&
    lbound(fortran_type_f%ts_intensity,2),&
    lbound(fortran_type_f%ts_intensity,3),&
    lbound(fortran_type_f%ts_intensity,4))
#ifdef ifort
  transfer_struct_c%ts_intensity_f_byte_size = transfer_struct_c%ts_intensity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_intensity_f_shapes(1) = size(fortran_type_f%ts_intensity, 1)
  transfer_struct_c%ts_intensity_f_shapes(2) = size(fortran_type_f%ts_intensity, 2)
  transfer_struct_c%ts_intensity_f_shapes(3) = size(fortran_type_f%ts_intensity, 3)
  transfer_struct_c%ts_intensity_f_shapes(4) = size(fortran_type_f%ts_intensity, 4)
  
  
  fortran_type_f%ts_mean_intensity = 0_fpk
  ts_mean_intensity_ptr => fortran_type_f%ts_mean_intensity
  transfer_struct_c%ts_mean_intensity = c_loc(ts_mean_intensity_ptr(&
    lbound(fortran_type_f%ts_mean_intensity,1),&
    lbound(fortran_type_f%ts_mean_intensity,2),&
    lbound(fortran_type_f%ts_mean_intensity,3),&
    lbound(fortran_type_f%ts_mean_intensity,4)))
  inquire(iolength=transfer_struct_c%ts_mean_intensity_f_byte_size) fortran_type_f%ts_mean_intensity(&
    lbound(fortran_type_f%ts_mean_intensity,1),&
    lbound(fortran_type_f%ts_mean_intensity,2),&
    lbound(fortran_type_f%ts_mean_intensity,3),&
    lbound(fortran_type_f%ts_mean_intensity,4))
#ifdef ifort
  transfer_struct_c%ts_mean_intensity_f_byte_size = transfer_struct_c%ts_mean_intensity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_mean_intensity_f_shapes(1) = size(fortran_type_f%ts_mean_intensity, 1)
  transfer_struct_c%ts_mean_intensity_f_shapes(2) = size(fortran_type_f%ts_mean_intensity, 2)
  transfer_struct_c%ts_mean_intensity_f_shapes(3) = size(fortran_type_f%ts_mean_intensity, 3)
  transfer_struct_c%ts_mean_intensity_f_shapes(4) = size(fortran_type_f%ts_mean_intensity, 4)
  
  
  fortran_type_f%ts_flux_integral = 0_fpk
  ts_flux_integral_ptr => fortran_type_f%ts_flux_integral
  transfer_struct_c%ts_flux_integral = c_loc(ts_flux_integral_ptr(&
    lbound(fortran_type_f%ts_flux_integral,1),&
    lbound(fortran_type_f%ts_flux_integral,2),&
    lbound(fortran_type_f%ts_flux_integral,3),&
    lbound(fortran_type_f%ts_flux_integral,4)))
  inquire(iolength=transfer_struct_c%ts_flux_integral_f_byte_size) fortran_type_f%ts_flux_integral(&
    lbound(fortran_type_f%ts_flux_integral,1),&
    lbound(fortran_type_f%ts_flux_integral,2),&
    lbound(fortran_type_f%ts_flux_integral,3),&
    lbound(fortran_type_f%ts_flux_integral,4))
#ifdef ifort
  transfer_struct_c%ts_flux_integral_f_byte_size = transfer_struct_c%ts_flux_integral_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_flux_integral_f_shapes(1) = size(fortran_type_f%ts_flux_integral, 1)
  transfer_struct_c%ts_flux_integral_f_shapes(2) = size(fortran_type_f%ts_flux_integral, 2)
  transfer_struct_c%ts_flux_integral_f_shapes(3) = size(fortran_type_f%ts_flux_integral, 3)
  transfer_struct_c%ts_flux_integral_f_shapes(4) = size(fortran_type_f%ts_flux_integral, 4)
  
  
  fortran_type_f%ts_dnflux_direct = 0_fpk
  ts_dnflux_direct_ptr => fortran_type_f%ts_dnflux_direct
  transfer_struct_c%ts_dnflux_direct = c_loc(ts_dnflux_direct_ptr(&
    lbound(fortran_type_f%ts_dnflux_direct,1),&
    lbound(fortran_type_f%ts_dnflux_direct,2),&
    lbound(fortran_type_f%ts_dnflux_direct,3)))
  inquire(iolength=transfer_struct_c%ts_dnflux_direct_f_byte_size) fortran_type_f%ts_dnflux_direct(&
    lbound(fortran_type_f%ts_dnflux_direct,1),&
    lbound(fortran_type_f%ts_dnflux_direct,2),&
    lbound(fortran_type_f%ts_dnflux_direct,3))
#ifdef ifort
  transfer_struct_c%ts_dnflux_direct_f_byte_size = transfer_struct_c%ts_dnflux_direct_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnflux_direct_f_shapes(1) = size(fortran_type_f%ts_dnflux_direct, 1)
  transfer_struct_c%ts_dnflux_direct_f_shapes(2) = size(fortran_type_f%ts_dnflux_direct, 2)
  transfer_struct_c%ts_dnflux_direct_f_shapes(3) = size(fortran_type_f%ts_dnflux_direct, 3)
  
  
  fortran_type_f%ts_dnmean_direct = 0_fpk
  ts_dnmean_direct_ptr => fortran_type_f%ts_dnmean_direct
  transfer_struct_c%ts_dnmean_direct = c_loc(ts_dnmean_direct_ptr(&
    lbound(fortran_type_f%ts_dnmean_direct,1),&
    lbound(fortran_type_f%ts_dnmean_direct,2),&
    lbound(fortran_type_f%ts_dnmean_direct,3)))
  inquire(iolength=transfer_struct_c%ts_dnmean_direct_f_byte_size) fortran_type_f%ts_dnmean_direct(&
    lbound(fortran_type_f%ts_dnmean_direct,1),&
    lbound(fortran_type_f%ts_dnmean_direct,2),&
    lbound(fortran_type_f%ts_dnmean_direct,3))
#ifdef ifort
  transfer_struct_c%ts_dnmean_direct_f_byte_size = transfer_struct_c%ts_dnmean_direct_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_dnmean_direct_f_shapes(1) = size(fortran_type_f%ts_dnmean_direct, 1)
  transfer_struct_c%ts_dnmean_direct_f_shapes(2) = size(fortran_type_f%ts_dnmean_direct, 2)
  transfer_struct_c%ts_dnmean_direct_f_shapes(3) = size(fortran_type_f%ts_dnmean_direct, 3)
  
  
  fortran_type_f%ts_fourier_saved = 0
  ts_fourier_saved_ptr => fortran_type_f%ts_fourier_saved
  transfer_struct_c%ts_fourier_saved = c_loc(ts_fourier_saved_ptr(&
    lbound(fortran_type_f%ts_fourier_saved,1),&
    lbound(fortran_type_f%ts_fourier_saved,2)))
  inquire(iolength=transfer_struct_c%ts_fourier_saved_f_byte_size) fortran_type_f%ts_fourier_saved(&
    lbound(fortran_type_f%ts_fourier_saved,1),&
    lbound(fortran_type_f%ts_fourier_saved,2))
#ifdef ifort
  transfer_struct_c%ts_fourier_saved_f_byte_size = transfer_struct_c%ts_fourier_saved_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_fourier_saved_f_shapes(1) = size(fortran_type_f%ts_fourier_saved, 1)
  transfer_struct_c%ts_fourier_saved_f_shapes(2) = size(fortran_type_f%ts_fourier_saved, 2)
  
  
  fortran_type_f%ts_n_geometries = 0
  ts_n_geometries_ptr => fortran_type_f%ts_n_geometries
  transfer_struct_c%ts_n_geometries = c_loc(ts_n_geometries_ptr)
  inquire(iolength=transfer_struct_c%ts_n_geometries_f_byte_size) fortran_type_f%ts_n_geometries
#ifdef ifort
  transfer_struct_c%ts_n_geometries_f_byte_size = transfer_struct_c%ts_n_geometries_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_main_outputs_c_init_only

subroutine lidort_main_outputs_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_main_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_main_outputs_c_destroy

subroutine lidort_main_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def, only : lidort_main_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_main_outputs), pointer :: fortran_type_f_from
  type(lidort_main_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_intensity = fortran_type_f_from%ts_intensity
  fortran_type_f_to%ts_mean_intensity = fortran_type_f_from%ts_mean_intensity
  fortran_type_f_to%ts_flux_integral = fortran_type_f_from%ts_flux_integral
  fortran_type_f_to%ts_dnflux_direct = fortran_type_f_from%ts_dnflux_direct
  fortran_type_f_to%ts_dnmean_direct = fortran_type_f_from%ts_dnmean_direct
  fortran_type_f_to%ts_fourier_saved = fortran_type_f_from%ts_fourier_saved
  fortran_type_f_to%ts_n_geometries = fortran_type_f_from%ts_n_geometries
  

end subroutine lidort_main_outputs_c_copy

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

  type(lidort_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_exception_handling_c_alloc_init

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def
  use lidort_pars

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
  use lidort_outputs_def, only : lidort_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_exception_handling_c_destroy

subroutine lidort_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_input_exception_handling_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

  type(lidort_input_exception_handling_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_input_exception_handling_c_alloc_init

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_input_exception_handling_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def
  use lidort_pars

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
  use lidort_outputs_def, only : lidort_input_exception_handling

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_input_exception_handling_c_destroy

subroutine lidort_input_exception_handling_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

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

! Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_outputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_outputs

  type(lidort_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_outputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_outputs_c_alloc_init

! Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Initializes only with no allocation
subroutine lidort_outputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_outputs_def
  use lidort_pars

  type(lidort_outputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  type(lidort_main_outputs), pointer :: main_ptr
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
  
  
  status_ptr => fortran_type_f%status
  transfer_struct_c%status = c_loc(status_ptr)
  inquire(iolength=transfer_struct_c%status_f_byte_size) fortran_type_f%status
#ifdef ifort
  transfer_struct_c%status_f_byte_size = transfer_struct_c%status_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_outputs_c_init_only

subroutine lidort_outputs_c_destroy(fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_outputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_outputs_c_destroy

subroutine lidort_outputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_outputs_def, only : lidort_outputs

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_outputs), pointer :: fortran_type_f_from
  type(lidort_outputs), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%main = fortran_type_f_from%main
  fortran_type_f_to%status = fortran_type_f_from%status
  

end subroutine lidort_outputs_c_copy

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_sup_brdf_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_brdf_def, only : lidort_sup_brdf

  type(lidort_sup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_brdf_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_brdf_c_alloc_init

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_brdf_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_brdf_def
  use lidort_pars

  type(lidort_sup_brdf_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_exactdb_brdfunc_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_brdf_f_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_user_brdf_f_0_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_user_brdf_f_ptr
  real(c_double), dimension(:,:), pointer :: ts_emissivity_ptr
  real(c_double), dimension(:,:), pointer :: ts_user_emissivity_ptr
  

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
    lbound(fortran_type_f%ts_emissivity,1),&
    lbound(fortran_type_f%ts_emissivity,2)))
  inquire(iolength=transfer_struct_c%ts_emissivity_f_byte_size) fortran_type_f%ts_emissivity(&
    lbound(fortran_type_f%ts_emissivity,1),&
    lbound(fortran_type_f%ts_emissivity,2))
#ifdef ifort
  transfer_struct_c%ts_emissivity_f_byte_size = transfer_struct_c%ts_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_emissivity_f_shapes(1) = size(fortran_type_f%ts_emissivity, 1)
  transfer_struct_c%ts_emissivity_f_shapes(2) = size(fortran_type_f%ts_emissivity, 2)
  
  
  fortran_type_f%ts_user_emissivity = 0_fpk
  ts_user_emissivity_ptr => fortran_type_f%ts_user_emissivity
  transfer_struct_c%ts_user_emissivity = c_loc(ts_user_emissivity_ptr(&
    lbound(fortran_type_f%ts_user_emissivity,1),&
    lbound(fortran_type_f%ts_user_emissivity,2)))
  inquire(iolength=transfer_struct_c%ts_user_emissivity_f_byte_size) fortran_type_f%ts_user_emissivity(&
    lbound(fortran_type_f%ts_user_emissivity,1),&
    lbound(fortran_type_f%ts_user_emissivity,2))
#ifdef ifort
  transfer_struct_c%ts_user_emissivity_f_byte_size = transfer_struct_c%ts_user_emissivity_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_user_emissivity_f_shapes(1) = size(fortran_type_f%ts_user_emissivity, 1)
  transfer_struct_c%ts_user_emissivity_f_shapes(2) = size(fortran_type_f%ts_user_emissivity, 2)
  
  
end subroutine lidort_sup_brdf_c_init_only

subroutine lidort_sup_brdf_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_brdf_def, only : lidort_sup_brdf

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_brdf_c_destroy

subroutine lidort_sup_brdf_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_brdf_def, only : lidort_sup_brdf

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

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_sup_sleave_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_sleave_def, only : lidort_sup_sleave

  type(lidort_sup_sleave_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_sleave_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_sleave_c_alloc_init

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_sleave_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_sleave_def
  use lidort_pars

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
  use lidort_sup_sleave_def, only : lidort_sup_sleave

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_sleave_c_destroy

subroutine lidort_sup_sleave_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_sleave_def, only : lidort_sup_sleave

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

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_sup_ss_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_ss_def, only : lidort_sup_ss

  type(lidort_sup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_ss_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_ss_c_alloc_init

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_ss_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_ss_def
  use lidort_pars

  type(lidort_sup_ss_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  real(c_double), dimension(:,:,:), pointer :: ts_intensity_ss_ptr
  real(c_double), dimension(:,:), pointer :: ts_intensity_db_ptr
  

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
  
  
end subroutine lidort_sup_ss_c_init_only

subroutine lidort_sup_ss_c_destroy(fortran_type_c) bind(C)
  use lidort_sup_ss_def, only : lidort_sup_ss

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_ss_c_destroy

subroutine lidort_sup_ss_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_ss_def, only : lidort_sup_ss

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_sup_ss), pointer :: fortran_type_f_from
  type(lidort_sup_ss), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_intensity_ss = fortran_type_f_from%ts_intensity_ss
  fortran_type_f_to%ts_intensity_db = fortran_type_f_from%ts_intensity_db
  

end subroutine lidort_sup_ss_c_copy

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
! Allocs and initializes type
subroutine lidort_sup_inout_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_inout_def, only : lidort_sup_inout

  type(lidort_sup_inout_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_sup_inout_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_sup_inout_c_alloc_init

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
! Initializes only with no allocation
subroutine lidort_sup_inout_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_sup_inout_def
  use lidort_sup_brdf_def
  use lidort_sup_ss_def
  use lidort_sup_sleave_def

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
  use lidort_sup_inout_def, only : lidort_sup_inout

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_sup_inout_c_destroy

subroutine lidort_sup_inout_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_sup_inout_def, only : lidort_sup_inout

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

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_boolean_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_boolean

  type(lidort_fixed_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_boolean_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_boolean_c_alloc_init

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_boolean_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_fixed_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_fullrad_mode_ptr
  logical(kind=4), pointer :: ts_do_sscorr_truncation_ptr
  logical(kind=4), pointer :: ts_do_ss_external_ptr
  logical(kind=4), pointer :: ts_do_ssfull_ptr
  logical(kind=4), pointer :: ts_do_thermal_emission_ptr
  logical(kind=4), pointer :: ts_do_surface_emission_ptr
  logical(kind=4), pointer :: ts_do_plane_parallel_ptr
  logical(kind=4), pointer :: ts_do_brdf_surface_ptr
  logical(kind=4), pointer :: ts_do_upwelling_ptr
  logical(kind=4), pointer :: ts_do_dnwelling_ptr
  logical(kind=4), pointer :: ts_do_surface_leaving_ptr
  logical(kind=4), pointer :: ts_do_sl_isotropic_ptr
  

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
  
  
  
  fortran_type_f%ts_do_sscorr_truncation = .FALSE.
  ts_do_sscorr_truncation_ptr => fortran_type_f%ts_do_sscorr_truncation
  transfer_struct_c%ts_do_sscorr_truncation = c_loc(ts_do_sscorr_truncation_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sscorr_truncation_f_byte_size) fortran_type_f%ts_do_sscorr_truncation
#ifdef ifort
  transfer_struct_c%ts_do_sscorr_truncation_f_byte_size = transfer_struct_c%ts_do_sscorr_truncation_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_ss_external = .FALSE.
  ts_do_ss_external_ptr => fortran_type_f%ts_do_ss_external
  transfer_struct_c%ts_do_ss_external = c_loc(ts_do_ss_external_ptr)
  inquire(iolength=transfer_struct_c%ts_do_ss_external_f_byte_size) fortran_type_f%ts_do_ss_external
#ifdef ifort
  transfer_struct_c%ts_do_ss_external_f_byte_size = transfer_struct_c%ts_do_ss_external_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_ssfull = .FALSE.
  ts_do_ssfull_ptr => fortran_type_f%ts_do_ssfull
  transfer_struct_c%ts_do_ssfull = c_loc(ts_do_ssfull_ptr)
  inquire(iolength=transfer_struct_c%ts_do_ssfull_f_byte_size) fortran_type_f%ts_do_ssfull
#ifdef ifort
  transfer_struct_c%ts_do_ssfull_f_byte_size = transfer_struct_c%ts_do_ssfull_f_byte_size * 4
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
  
  
  
end subroutine lidort_fixed_boolean_c_init_only

subroutine lidort_fixed_boolean_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_boolean

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_boolean_c_destroy

subroutine lidort_fixed_boolean_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_boolean

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_boolean), pointer :: fortran_type_f_from
  type(lidort_fixed_boolean), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_fullrad_mode = fortran_type_f_from%ts_do_fullrad_mode
  fortran_type_f_to%ts_do_sscorr_truncation = fortran_type_f_from%ts_do_sscorr_truncation
  fortran_type_f_to%ts_do_ss_external = fortran_type_f_from%ts_do_ss_external
  fortran_type_f_to%ts_do_ssfull = fortran_type_f_from%ts_do_ssfull
  fortran_type_f_to%ts_do_thermal_emission = fortran_type_f_from%ts_do_thermal_emission
  fortran_type_f_to%ts_do_surface_emission = fortran_type_f_from%ts_do_surface_emission
  fortran_type_f_to%ts_do_plane_parallel = fortran_type_f_from%ts_do_plane_parallel
  fortran_type_f_to%ts_do_brdf_surface = fortran_type_f_from%ts_do_brdf_surface
  fortran_type_f_to%ts_do_upwelling = fortran_type_f_from%ts_do_upwelling
  fortran_type_f_to%ts_do_dnwelling = fortran_type_f_from%ts_do_dnwelling
  fortran_type_f_to%ts_do_surface_leaving = fortran_type_f_from%ts_do_surface_leaving
  fortran_type_f_to%ts_do_sl_isotropic = fortran_type_f_from%ts_do_sl_isotropic
  

end subroutine lidort_fixed_boolean_c_copy

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_control_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_control

  type(lidort_fixed_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_control_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_control_c_alloc_init

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_control_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_fixed_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_nstreams_ptr
  integer(c_int), pointer :: ts_nlayers_ptr
  integer(c_int), pointer :: ts_nfinelayers_ptr
  integer(c_int), pointer :: ts_n_thermal_coeffs_ptr
  real(c_double), pointer :: ts_lidort_accuracy_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
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
  
  
  
end subroutine lidort_fixed_control_c_init_only

subroutine lidort_fixed_control_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_control

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_control_c_destroy

subroutine lidort_fixed_control_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_control

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_control), pointer :: fortran_type_f_from
  type(lidort_fixed_control), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_nstreams = fortran_type_f_from%ts_nstreams
  fortran_type_f_to%ts_nlayers = fortran_type_f_from%ts_nlayers
  fortran_type_f_to%ts_nfinelayers = fortran_type_f_from%ts_nfinelayers
  fortran_type_f_to%ts_n_thermal_coeffs = fortran_type_f_from%ts_n_thermal_coeffs
  fortran_type_f_to%ts_lidort_accuracy = fortran_type_f_from%ts_lidort_accuracy
  

end subroutine lidort_fixed_control_c_copy

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_sunrays_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_sunrays

  type(lidort_fixed_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_sunrays_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_sunrays_c_alloc_init

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_sunrays_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_fixed_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_sunrays_c_destroy

subroutine lidort_fixed_sunrays_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_sunrays), pointer :: fortran_type_f_from
  type(lidort_fixed_sunrays), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_flux_factor = fortran_type_f_from%ts_flux_factor
  

end subroutine lidort_fixed_sunrays_c_copy

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_uservalues_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_uservalues

  type(lidort_fixed_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_uservalues_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_uservalues_c_alloc_init

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_uservalues_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_fixed_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_n_user_streams_ptr
  integer(c_int), pointer :: ts_n_user_levels_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_n_user_streams = 0
  ts_n_user_streams_ptr => fortran_type_f%ts_n_user_streams
  transfer_struct_c%ts_n_user_streams = c_loc(ts_n_user_streams_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_streams_f_byte_size) fortran_type_f%ts_n_user_streams
#ifdef ifort
  transfer_struct_c%ts_n_user_streams_f_byte_size = transfer_struct_c%ts_n_user_streams_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_n_user_levels = 0
  ts_n_user_levels_ptr => fortran_type_f%ts_n_user_levels
  transfer_struct_c%ts_n_user_levels = c_loc(ts_n_user_levels_ptr)
  inquire(iolength=transfer_struct_c%ts_n_user_levels_f_byte_size) fortran_type_f%ts_n_user_levels
#ifdef ifort
  transfer_struct_c%ts_n_user_levels_f_byte_size = transfer_struct_c%ts_n_user_levels_f_byte_size * 4
#endif
  
  
  
end subroutine lidort_fixed_uservalues_c_init_only

subroutine lidort_fixed_uservalues_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_uservalues_c_destroy

subroutine lidort_fixed_uservalues_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_uservalues), pointer :: fortran_type_f_from
  type(lidort_fixed_uservalues), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_n_user_streams = fortran_type_f_from%ts_n_user_streams
  fortran_type_f_to%ts_n_user_levels = fortran_type_f_from%ts_n_user_levels
  

end subroutine lidort_fixed_uservalues_c_copy

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_chapman_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_chapman

  type(lidort_fixed_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_chapman_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_chapman_c_alloc_init

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_chapman_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_fixed_chapman

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_chapman_c_destroy

subroutine lidort_fixed_chapman_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_chapman

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

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_optical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_optical

  type(lidort_fixed_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_optical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_optical_c_alloc_init

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_optical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_fixed_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  real(c_double), dimension(:,:), pointer :: ts_deltau_vert_input_ptr
  real(c_double), dimension(:,:,:), pointer :: ts_phasmoms_total_input_ptr
  real(c_double), dimension(:,:), pointer :: ts_thermal_bb_input_ptr
  real(c_double), dimension(:), pointer :: ts_lambertian_albedo_ptr
  real(c_double), dimension(:), pointer :: ts_surface_bb_input_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_deltau_vert_input = 0_fpk
  ts_deltau_vert_input_ptr => fortran_type_f%ts_deltau_vert_input
  transfer_struct_c%ts_deltau_vert_input = c_loc(ts_deltau_vert_input_ptr(&
    lbound(fortran_type_f%ts_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_deltau_vert_input,2)))
  inquire(iolength=transfer_struct_c%ts_deltau_vert_input_f_byte_size) fortran_type_f%ts_deltau_vert_input(&
    lbound(fortran_type_f%ts_deltau_vert_input,1),&
    lbound(fortran_type_f%ts_deltau_vert_input,2))
#ifdef ifort
  transfer_struct_c%ts_deltau_vert_input_f_byte_size = transfer_struct_c%ts_deltau_vert_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_deltau_vert_input_f_shapes(1) = size(fortran_type_f%ts_deltau_vert_input, 1)
  transfer_struct_c%ts_deltau_vert_input_f_shapes(2) = size(fortran_type_f%ts_deltau_vert_input, 2)
  
  
  fortran_type_f%ts_phasmoms_total_input = 0_fpk
  ts_phasmoms_total_input_ptr => fortran_type_f%ts_phasmoms_total_input
  transfer_struct_c%ts_phasmoms_total_input = c_loc(ts_phasmoms_total_input_ptr(&
    lbound(fortran_type_f%ts_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_phasmoms_total_input,3)))
  inquire(iolength=transfer_struct_c%ts_phasmoms_total_input_f_byte_size) fortran_type_f%ts_phasmoms_total_input(&
    lbound(fortran_type_f%ts_phasmoms_total_input,1),&
    lbound(fortran_type_f%ts_phasmoms_total_input,2),&
    lbound(fortran_type_f%ts_phasmoms_total_input,3))
#ifdef ifort
  transfer_struct_c%ts_phasmoms_total_input_f_byte_size = transfer_struct_c%ts_phasmoms_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_phasmoms_total_input_f_shapes(1) = size(fortran_type_f%ts_phasmoms_total_input, 1)
  transfer_struct_c%ts_phasmoms_total_input_f_shapes(2) = size(fortran_type_f%ts_phasmoms_total_input, 2)
  transfer_struct_c%ts_phasmoms_total_input_f_shapes(3) = size(fortran_type_f%ts_phasmoms_total_input, 3)
  
  
  fortran_type_f%ts_thermal_bb_input = 0_fpk
  ts_thermal_bb_input_ptr => fortran_type_f%ts_thermal_bb_input
  transfer_struct_c%ts_thermal_bb_input = c_loc(ts_thermal_bb_input_ptr(&
    lbound(fortran_type_f%ts_thermal_bb_input,1),&
    lbound(fortran_type_f%ts_thermal_bb_input,2)))
  inquire(iolength=transfer_struct_c%ts_thermal_bb_input_f_byte_size) fortran_type_f%ts_thermal_bb_input(&
    lbound(fortran_type_f%ts_thermal_bb_input,1),&
    lbound(fortran_type_f%ts_thermal_bb_input,2))
#ifdef ifort
  transfer_struct_c%ts_thermal_bb_input_f_byte_size = transfer_struct_c%ts_thermal_bb_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_thermal_bb_input_f_shapes(1) = size(fortran_type_f%ts_thermal_bb_input, 1)
  transfer_struct_c%ts_thermal_bb_input_f_shapes(2) = size(fortran_type_f%ts_thermal_bb_input, 2)
  
  
  fortran_type_f%ts_lambertian_albedo = 0_fpk
  ts_lambertian_albedo_ptr => fortran_type_f%ts_lambertian_albedo
  transfer_struct_c%ts_lambertian_albedo = c_loc(ts_lambertian_albedo_ptr(&
    lbound(fortran_type_f%ts_lambertian_albedo,1)))
  inquire(iolength=transfer_struct_c%ts_lambertian_albedo_f_byte_size) fortran_type_f%ts_lambertian_albedo(&
    lbound(fortran_type_f%ts_lambertian_albedo,1))
#ifdef ifort
  transfer_struct_c%ts_lambertian_albedo_f_byte_size = transfer_struct_c%ts_lambertian_albedo_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_lambertian_albedo_f_shapes(1) = size(fortran_type_f%ts_lambertian_albedo, 1)
  
  
  fortran_type_f%ts_surface_bb_input = 0_fpk
  ts_surface_bb_input_ptr => fortran_type_f%ts_surface_bb_input
  transfer_struct_c%ts_surface_bb_input = c_loc(ts_surface_bb_input_ptr(&
    lbound(fortran_type_f%ts_surface_bb_input,1)))
  inquire(iolength=transfer_struct_c%ts_surface_bb_input_f_byte_size) fortran_type_f%ts_surface_bb_input(&
    lbound(fortran_type_f%ts_surface_bb_input,1))
#ifdef ifort
  transfer_struct_c%ts_surface_bb_input_f_byte_size = transfer_struct_c%ts_surface_bb_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_surface_bb_input_f_shapes(1) = size(fortran_type_f%ts_surface_bb_input, 1)
  
  
end subroutine lidort_fixed_optical_c_init_only

subroutine lidort_fixed_optical_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_optical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_optical_c_destroy

subroutine lidort_fixed_optical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_optical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_fixed_optical), pointer :: fortran_type_f_from
  type(lidort_fixed_optical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_deltau_vert_input = fortran_type_f_from%ts_deltau_vert_input
  fortran_type_f_to%ts_phasmoms_total_input = fortran_type_f_from%ts_phasmoms_total_input
  fortran_type_f_to%ts_thermal_bb_input = fortran_type_f_from%ts_thermal_bb_input
  fortran_type_f_to%ts_lambertian_albedo = fortran_type_f_from%ts_lambertian_albedo
  fortran_type_f_to%ts_surface_bb_input = fortran_type_f_from%ts_surface_bb_input
  

end subroutine lidort_fixed_optical_c_copy

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_inputs

  type(lidort_fixed_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_fixed_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_fixed_inputs_c_alloc_init

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_fixed_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_fixed_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  type(lidort_fixed_boolean), pointer :: bool_ptr
  type(lidort_fixed_control), pointer :: cont_ptr
  type(lidort_fixed_sunrays), pointer :: sunrays_ptr
  type(lidort_fixed_uservalues), pointer :: userval_ptr
  type(lidort_fixed_chapman), pointer :: chapman_ptr
  type(lidort_fixed_optical), pointer :: optical_ptr
  

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
  
  
  
end subroutine lidort_fixed_inputs_c_init_only

subroutine lidort_fixed_inputs_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_fixed_inputs_c_destroy

subroutine lidort_fixed_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_fixed_inputs

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
  

end subroutine lidort_fixed_inputs_c_copy

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_boolean_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_boolean

  type(lidort_modified_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_boolean_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_boolean_c_alloc_init

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_boolean_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_modified_boolean_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  logical(kind=4), pointer :: ts_do_sscorr_nadir_ptr
  logical(kind=4), pointer :: ts_do_sscorr_outgoing_ptr
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
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_do_sscorr_nadir = .FALSE.
  ts_do_sscorr_nadir_ptr => fortran_type_f%ts_do_sscorr_nadir
  transfer_struct_c%ts_do_sscorr_nadir = c_loc(ts_do_sscorr_nadir_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sscorr_nadir_f_byte_size) fortran_type_f%ts_do_sscorr_nadir
#ifdef ifort
  transfer_struct_c%ts_do_sscorr_nadir_f_byte_size = transfer_struct_c%ts_do_sscorr_nadir_f_byte_size * 4
#endif
  
  
  
  fortran_type_f%ts_do_sscorr_outgoing = .FALSE.
  ts_do_sscorr_outgoing_ptr => fortran_type_f%ts_do_sscorr_outgoing
  transfer_struct_c%ts_do_sscorr_outgoing = c_loc(ts_do_sscorr_outgoing_ptr)
  inquire(iolength=transfer_struct_c%ts_do_sscorr_outgoing_f_byte_size) fortran_type_f%ts_do_sscorr_outgoing
#ifdef ifort
  transfer_struct_c%ts_do_sscorr_outgoing_f_byte_size = transfer_struct_c%ts_do_sscorr_outgoing_f_byte_size * 4
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
  
  
  
end subroutine lidort_modified_boolean_c_init_only

subroutine lidort_modified_boolean_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_boolean

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_boolean_c_destroy

subroutine lidort_modified_boolean_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_boolean

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_boolean), pointer :: fortran_type_f_from
  type(lidort_modified_boolean), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_do_sscorr_nadir = fortran_type_f_from%ts_do_sscorr_nadir
  fortran_type_f_to%ts_do_sscorr_outgoing = fortran_type_f_from%ts_do_sscorr_outgoing
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
  

end subroutine lidort_modified_boolean_c_copy

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_control_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_control

  type(lidort_modified_control_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_control_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_control_c_alloc_init

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_control_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_modified_control

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_control_c_destroy

subroutine lidort_modified_control_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_control

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_control), pointer :: fortran_type_f_from
  type(lidort_modified_control), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_nmoments_input = fortran_type_f_from%ts_nmoments_input
  

end subroutine lidort_modified_control_c_copy

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_sunrays_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_sunrays

  type(lidort_modified_sunrays_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_sunrays_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_sunrays_c_alloc_init

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_sunrays_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_modified_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_sunrays_c_destroy

subroutine lidort_modified_sunrays_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_sunrays

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_sunrays), pointer :: fortran_type_f_from
  type(lidort_modified_sunrays), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_nbeams = fortran_type_f_from%ts_nbeams
  fortran_type_f_to%ts_beam_szas = fortran_type_f_from%ts_beam_szas
  

end subroutine lidort_modified_sunrays_c_copy

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_uservalues_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_uservalues

  type(lidort_modified_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_uservalues_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_uservalues_c_alloc_init

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_uservalues_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_modified_uservalues_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  integer(c_int), pointer :: ts_n_user_relazms_ptr
  real(c_double), dimension(:), pointer :: ts_user_relazms_ptr
  real(c_double), dimension(:), pointer :: ts_user_angles_input_ptr
  real(c_double), dimension(:), pointer :: ts_user_levels_ptr
  real(c_double), pointer :: ts_geometry_specheight_ptr
  

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
  
  
  
end subroutine lidort_modified_uservalues_c_init_only

subroutine lidort_modified_uservalues_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_uservalues_c_destroy

subroutine lidort_modified_uservalues_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_uservalues

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_uservalues), pointer :: fortran_type_f_from
  type(lidort_modified_uservalues), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_n_user_relazms = fortran_type_f_from%ts_n_user_relazms
  fortran_type_f_to%ts_user_relazms = fortran_type_f_from%ts_user_relazms
  fortran_type_f_to%ts_user_angles_input = fortran_type_f_from%ts_user_angles_input
  fortran_type_f_to%ts_user_levels = fortran_type_f_from%ts_user_levels
  fortran_type_f_to%ts_geometry_specheight = fortran_type_f_from%ts_geometry_specheight
  

end subroutine lidort_modified_uservalues_c_copy

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_chapman_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_chapman

  type(lidort_modified_chapman_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_chapman_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_chapman_c_alloc_init

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_chapman_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_modified_chapman

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_chapman_c_destroy

subroutine lidort_modified_chapman_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_chapman

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_chapman), pointer :: fortran_type_f_from
  type(lidort_modified_chapman), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_earth_radius = fortran_type_f_from%ts_earth_radius
  

end subroutine lidort_modified_chapman_c_copy

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_optical_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_optical

  type(lidort_modified_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_optical_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_optical_c_alloc_init

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_optical_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

  type(lidort_modified_optical_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  real(c_double), dimension(:,:), pointer :: ts_omega_total_input_ptr
  

  call c_f_pointer(fortran_type_c, fortran_type_f)

  ! For each variable in the structure:
  !   * Initalize
  !   * Set pointer location into C transfer struct
  !   * Set bit size and extents

  
  fortran_type_f%ts_omega_total_input = 0_fpk
  ts_omega_total_input_ptr => fortran_type_f%ts_omega_total_input
  transfer_struct_c%ts_omega_total_input = c_loc(ts_omega_total_input_ptr(&
    lbound(fortran_type_f%ts_omega_total_input,1),&
    lbound(fortran_type_f%ts_omega_total_input,2)))
  inquire(iolength=transfer_struct_c%ts_omega_total_input_f_byte_size) fortran_type_f%ts_omega_total_input(&
    lbound(fortran_type_f%ts_omega_total_input,1),&
    lbound(fortran_type_f%ts_omega_total_input,2))
#ifdef ifort
  transfer_struct_c%ts_omega_total_input_f_byte_size = transfer_struct_c%ts_omega_total_input_f_byte_size * 4
#endif
  
  transfer_struct_c%ts_omega_total_input_f_shapes(1) = size(fortran_type_f%ts_omega_total_input, 1)
  transfer_struct_c%ts_omega_total_input_f_shapes(2) = size(fortran_type_f%ts_omega_total_input, 2)
  
  
end subroutine lidort_modified_optical_c_init_only

subroutine lidort_modified_optical_c_destroy(fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_optical

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_optical_c_destroy

subroutine lidort_modified_optical_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_optical

  type(c_ptr), intent(inout) :: fortran_type_c_from
  type(c_ptr), intent(inout) :: fortran_type_c_to

  type(lidort_modified_optical), pointer :: fortran_type_f_from
  type(lidort_modified_optical), pointer :: fortran_type_f_to

  call c_f_pointer(fortran_type_c_from, fortran_type_f_from)
  call c_f_pointer(fortran_type_c_to, fortran_type_f_to)

  fortran_type_f_to%ts_omega_total_input = fortran_type_f_from%ts_omega_total_input
  

end subroutine lidort_modified_optical_c_copy

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_inputs_c_alloc_init(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_inputs

  type(lidort_modified_inputs_c), intent(inout) :: transfer_struct_c
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  allocate(fortran_type_f)
  fortran_type_c = c_loc(fortran_type_f)

  call lidort_modified_inputs_c_init_only(transfer_struct_c, fortran_type_c) 

end subroutine lidort_modified_inputs_c_alloc_init

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Initializes only with no allocation
subroutine lidort_modified_inputs_c_init_only(transfer_struct_c, fortran_type_c) bind(C)
  use lidort_inputs_def
  use lidort_pars

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
  use lidort_inputs_def, only : lidort_modified_inputs

  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  deallocate(fortran_type_f)
end subroutine lidort_modified_inputs_c_destroy

subroutine lidort_modified_inputs_c_copy(fortran_type_c_from, fortran_type_c_to) bind(C)
  use lidort_inputs_def, only : lidort_modified_inputs

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


! Wrapper for character variable "bs_brdf_names" of type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
subroutine brdf_sup_inputs_bs_brdf_names_get(fortran_type_c, bs_brdf_names_in_shape_1, &
      bs_brdf_names_in_len, &
      bs_brdf_names_in) bind(C)
  use brdf_sup_inputs_def, only : brdf_sup_inputs

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
! Wrapper for character variable "bs_inputmessages" of type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
subroutine brdf_input_exception_handling_bs_inputmessages_get(fortran_type_c, bs_inputmessages_in_shape_1, &
      bs_inputmessages_in_len, &
      bs_inputmessages_in) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

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
! Wrapper for character variable "bs_inputactions" of type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
subroutine brdf_input_exception_handling_bs_inputactions_get(fortran_type_c, bs_inputactions_in_shape_1, &
      bs_inputactions_in_len, &
      bs_inputactions_in) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

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
! Wrapper for character variable "ts_checkmessages" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_checkmessages_get(fortran_type_c, ts_checkmessages_in_shape_1, &
      ts_checkmessages_in_len, &
      ts_checkmessages_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_actions" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_actions_get(fortran_type_c, ts_actions_in_shape_1, &
      ts_actions_in_len, &
      ts_actions_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_message" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_message_get(fortran_type_c, ts_message_in_len, &
      ts_message_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_trace_1" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_1_get(fortran_type_c, ts_trace_1_in_len, &
      ts_trace_1_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_trace_2" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_2_get(fortran_type_c, ts_trace_2_in_len, &
      ts_trace_2_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_trace_3" of type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine exception_handling_ts_trace_3_get(fortran_type_c, ts_trace_3_in_len, &
      ts_trace_3_in) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

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
! Wrapper for character variable "ts_inputmessages" of type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine input_exception_handling_ts_inputmessages_get(fortran_type_c, ts_inputmessages_in_shape_1, &
      ts_inputmessages_in_len, &
      ts_inputmessages_in) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

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
! Wrapper for character variable "ts_inputactions" of type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
subroutine input_exception_handling_ts_inputactions_get(fortran_type_c, ts_inputactions_in_shape_1, &
      ts_inputactions_in_len, &
      ts_inputactions_in) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

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


end module lidort_interface_types