
module TWOSTREAM_LS_BRDF_SUPPLEMENT_M_WRAP

use iso_c_binding
use twostream_ls_brdf_supplement_m
use twostream_brdf_supplement_m
use twostream_ls_brdfkernels_m

! This module was auto-generated 

implicit none

contains

subroutine twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap (maxbeams, &
                                                                        max_user_streams, &
                                                                        max_user_obsgeoms, &
                                                                        maxstreams_brdf, &
                                                                        max_brdf_kernels, &
                                                                        max_brdf_parameters, &
                                                                        max_surfacewfs, &
                                                                        do_solar_sources, &
                                                                        do_user_obsgeoms, &
                                                                        lambertian_kernel_flag, &
                                                                        do_shadow_effect, &
                                                                        do_surface_emission, &
                                                                        nbeams, &
                                                                        n_user_streams, &
                                                                        n_user_obsgeoms, &
                                                                        beam_szas, &
                                                                        user_angles, &
                                                                        user_obsgeoms, &
                                                                        stream_value, &
                                                                        nstreams_brdf, &
                                                                        n_brdf_kernels, &
                                                                        which_brdf, &
                                                                        brdf_factors, &
                                                                        n_brdf_parameters, &
                                                                        brdf_parameters, &
                                                                        do_kernel_factor_wfs, &
                                                                        do_kernel_params_wfs, &
                                                                        do_kparams_derivs, &
                                                                        n_surface_wfs, &
                                                                        n_kernel_factor_wfs, &
                                                                        n_kernel_params_wfs, &
                                                                        brdf_f_0, &
                                                                        brdf_f, &
                                                                        ubrdf_f, &
                                                                        emissivity, &
                                                                        ls_brdf_f_0, &
                                                                        ls_brdf_f, &
                                                                        ls_ubrdf_f, &
                                                                        ls_emissivity, &
                                                                        status_brdfsup, &
                                                                        message_len, &
                                                                        message, &
                                                                        action_len, &
                                                                        action) bind(C)
  use twostream_brdf_supplement_m
  use twostream_ls_brdfkernels_m

  ! Arguments
  integer(c_int), intent(in) :: maxbeams
  integer(c_int), intent(in) :: max_user_streams
  integer(c_int), intent(in) :: max_user_obsgeoms
  integer(c_int), intent(in) :: maxstreams_brdf
  integer(c_int), intent(in) :: max_brdf_kernels
  integer(c_int), intent(in) :: max_brdf_parameters
  integer(c_int), intent(in) :: max_surfacewfs
  logical(c_bool), intent(in) :: do_solar_sources
  logical(c_bool), intent(in) :: do_user_obsgeoms
  logical(c_bool), dimension(MAX_BRDF_KERNELS), intent(in) :: lambertian_kernel_flag
  logical(c_bool), intent(in) :: do_shadow_effect
  logical(c_bool), intent(in) :: do_surface_emission
  integer(c_int), intent(inout) :: nbeams
  integer(c_int), intent(inout) :: n_user_streams
  integer(c_int), intent(in) :: n_user_obsgeoms
  real(c_double), dimension(MAXBEAMS), intent(inout) :: beam_szas
  real(c_double), dimension(MAX_USER_STREAMS), intent(inout) :: user_angles
  real(c_double), dimension(MAX_USER_OBSGEOMS, 3), intent(in) :: user_obsgeoms
  real(c_double), intent(in) :: stream_value
  integer(c_int), intent(in) :: nstreams_brdf
  integer(c_int), intent(in) :: n_brdf_kernels
  integer(c_int), dimension(MAX_BRDF_KERNELS), intent(in) :: which_brdf
  real(c_double), dimension(MAX_BRDF_KERNELS), intent(in) :: brdf_factors
  integer(c_int), dimension(MAX_BRDF_KERNELS), intent(in) :: n_brdf_parameters
  real(c_double), dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS), intent(in) :: brdf_parameters
  logical(c_bool), dimension(MAX_BRDF_KERNELS), intent(in) :: do_kernel_factor_wfs
  logical(c_bool), dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS), intent(in) :: do_kernel_params_wfs
  logical(c_bool), dimension(MAX_BRDF_KERNELS), intent(out) :: do_kparams_derivs
  integer(c_int), intent(out) :: n_surface_wfs
  integer(c_int), intent(out) :: n_kernel_factor_wfs
  integer(c_int), intent(out) :: n_kernel_params_wfs
  real(c_double), dimension(0:1, MAXBEAMS), intent(out) :: brdf_f_0
  real(c_double), dimension(0:1), intent(out) :: brdf_f
  real(c_double), dimension(0:1, MAX_USER_STREAMS), intent(out) :: ubrdf_f
  real(c_double), intent(out) :: emissivity
  real(c_double), dimension(MAX_SURFACEWFS, 0:1, MAXBEAMS), intent(out) :: ls_brdf_f_0
  real(c_double), dimension(MAX_SURFACEWFS, 0:1), intent(out) :: ls_brdf_f
  real(c_double), dimension(MAX_SURFACEWFS, 0:1, MAX_USER_STREAMS), intent(out) :: ls_ubrdf_f
  real(c_double), dimension(MAX_SURFACEWFS), intent(out) :: ls_emissivity
  integer(c_int), intent(out) :: status_brdfsup
  integer(c_int), intent(in) :: message_len
  character(kind=c_char) , intent(inout) :: message(message_len+1)
  integer(c_int), intent(in) :: action_len
  character(kind=c_char) , intent(inout) :: action(action_len+1)

  ! Local variables
  logical(kind=4) :: do_solar_sources_lcl
  logical(kind=4) :: do_user_obsgeoms_lcl
  logical(kind=4), dimension(MAX_BRDF_KERNELS) :: lambertian_kernel_flag_lcl
  logical(kind=4) :: do_shadow_effect_lcl
  logical(kind=4) :: do_surface_emission_lcl
  logical(kind=4), dimension(MAX_BRDF_KERNELS) :: do_kernel_factor_wfs_lcl
  logical(kind=4), dimension(MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS) :: do_kernel_params_wfs_lcl
  logical(kind=4), dimension(MAX_BRDF_KERNELS) :: do_kparams_derivs_lcl
  character(kind=c_char, len=message_len) :: message_lcl
  integer :: len_idx
  character(kind=c_char, len=action_len) :: action_lcl

  ! Convert input arguments
  do_solar_sources_lcl = do_solar_sources
  do_user_obsgeoms_lcl = do_user_obsgeoms
  lambertian_kernel_flag_lcl = lambertian_kernel_flag
  do_shadow_effect_lcl = do_shadow_effect
  do_surface_emission_lcl = do_surface_emission
  do_kernel_factor_wfs_lcl = do_kernel_factor_wfs
  do_kernel_params_wfs_lcl = do_kernel_params_wfs

  call twostream_ls_brdfmaster(maxbeams, &
                               max_user_streams, &
                               max_user_obsgeoms, &
                               maxstreams_brdf, &
                               max_brdf_kernels, &
                               max_brdf_parameters, &
                               max_surfacewfs, &
                               do_solar_sources_lcl, &
                               do_user_obsgeoms_lcl, &
                               lambertian_kernel_flag_lcl, &
                               do_shadow_effect_lcl, &
                               do_surface_emission_lcl, &
                               nbeams, &
                               n_user_streams, &
                               n_user_obsgeoms, &
                               beam_szas, &
                               user_angles, &
                               user_obsgeoms, &
                               stream_value, &
                               nstreams_brdf, &
                               n_brdf_kernels, &
                               which_brdf, &
                               brdf_factors, &
                               n_brdf_parameters, &
                               brdf_parameters, &
                               do_kernel_factor_wfs_lcl, &
                               do_kernel_params_wfs_lcl, &
                               do_kparams_derivs_lcl, &
                               n_surface_wfs, &
                               n_kernel_factor_wfs, &
                               n_kernel_params_wfs, &
                               brdf_f_0, &
                               brdf_f, &
                               ubrdf_f, &
                               emissivity, &
                               ls_brdf_f_0, &
                               ls_brdf_f, &
                               ls_ubrdf_f, &
                               ls_emissivity, &
                               status_brdfsup, &
                               message_lcl, &
                               action_lcl)

  ! Convert output arguments
  do_kparams_derivs = do_kparams_derivs_lcl
  do len_idx = 1, message_len
    message(len_idx) = &
      message_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(message_lcl(:))+1
  message(len_idx) = c_null_char
  do len_idx = 1, action_len
    action(len_idx) = &
      action_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(action_lcl(:))+1
  action(len_idx) = c_null_char

end subroutine twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap

end module TWOSTREAM_LS_BRDF_SUPPLEMENT_M_WRAP

module TWOSTREAM_LPS_MASTER_M_WRAP

use iso_c_binding
use twostream_lps_master_m
use twostream_taylor_m
use twostream_inputs_m
use twostream_writemodules_m
use twostream_miscsetups_m
use twostream_solutions_m
use twostream_bvproblem_m
use twostream_intensity_m
use twostream_thermalsup_m
use twostream_l_inputs_m
use twostream_l_writemodules_m
use twostream_lp_miscsetups_m
use twostream_la_solutions_m
use twostream_lp_bvproblem_m
use twostream_ls_bvproblem_m
use twostream_lp_jacobians_m
use twostream_ls_jacobians_m
use twostream_thermalsup_plus_m
use twostream_lssl_jacobians_m

! This module was auto-generated 

implicit none

contains

subroutine twostream_lps_master_m_twostream_lps_master_wrap (maxlayers, &
                                                             maxtotal, &
                                                             maxmessages, &
                                                             maxbeams, &
                                                             max_geometries, &
                                                             max_user_streams, &
                                                             max_user_relazms, &
                                                             max_user_obsgeoms, &
                                                             max_atmoswfs, &
                                                             max_surfacewfs, &
                                                             max_sleavewfs, &
                                                             do_upwelling, &
                                                             do_dnwelling, &
                                                             do_plane_parallel, &
                                                             do_2s_levelout, &
                                                             do_mvout_only, &
                                                             do_additional_mvout, &
                                                             do_solar_sources, &
                                                             do_thermal_emission, &
                                                             do_surface_emission, &
                                                             do_d2s_scaling, &
                                                             do_brdf_surface, &
                                                             do_user_obsgeoms, &
                                                             do_surface_leaving, &
                                                             do_sl_isotropic, &
                                                             do_pentadiag_inverse, &
                                                             bvpindex, &
                                                             bvpscalefactor, &
                                                             taylor_order, &
                                                             taylor_small, &
                                                             tcutoff, &
                                                             nlayers, &
                                                             ntotal, &
                                                             stream_value, &
                                                             n_user_obsgeoms, &
                                                             user_obsgeoms, &
                                                             n_user_streams, &
                                                             user_angles, &
                                                             n_user_relazms, &
                                                             user_relazms, &
                                                             flux_factor, &
                                                             nbeams, &
                                                             beam_szas, &
                                                             earth_radius, &
                                                             height_grid, &
                                                             deltau_input, &
                                                             omega_input, &
                                                             asymm_input, &
                                                             d2s_scaling, &
                                                             thermal_bb_input, &
                                                             lambertian_albedo, &
                                                             brdf_f_0, &
                                                             brdf_f, &
                                                             ubrdf_f, &
                                                             emissivity, &
                                                             surfbb, &
                                                             slterm_isotropic, &
                                                             slterm_f_0, &
                                                             do_profile_wfs, &
                                                             do_surface_wfs, &
                                                             do_sleave_wfs, &
                                                             layer_vary_flag, &
                                                             layer_vary_number, &
                                                             n_surface_wfs, &
                                                             n_sleave_wfs, &
                                                             lssl_slterm_isotropic, &
                                                             lssl_slterm_f_0, &
                                                             l_deltau_input, &
                                                             l_omega_input, &
                                                             l_asymm_input, &
                                                             l_d2s_scaling, &
                                                             ls_brdf_f_0, &
                                                             ls_brdf_f, &
                                                             ls_ubrdf_f, &
                                                             ls_emissivity, &
                                                             intensity_toa, &
                                                             profilewf_toa, &
                                                             surfacewf_toa, &
                                                             intensity_boa, &
                                                             profilewf_boa, &
                                                             surfacewf_boa, &
                                                             radlevel_up, &
                                                             radlevel_dn, &
                                                             n_geometries, &
                                                             profjaclevel_up, &
                                                             profjaclevel_dn, &
                                                             surfjaclevel_up, &
                                                             surfjaclevel_dn, &
                                                             fluxes_toa, &
                                                             profjacfluxes_toa, &
                                                             surfjacfluxes_toa, &
                                                             fluxes_boa, &
                                                             profjacfluxes_boa, &
                                                             surfjacfluxes_boa, &
                                                             status_inputcheck, &
                                                             c_nmessages, &
                                                             c_messages_shape_1, &
                                                             c_messages_len, &
                                                             c_messages, &
                                                             c_actions_shape_1, &
                                                             c_actions_len, &
                                                             c_actions, &
                                                             status_execution, &
                                                             e_message_len, &
                                                             e_message, &
                                                             e_trace_1_len, &
                                                             e_trace_1, &
                                                             e_trace_2_len, &
                                                             e_trace_2) bind(C)
  use twostream_taylor_m
  use twostream_inputs_m
  use twostream_writemodules_m
  use twostream_miscsetups_m
  use twostream_solutions_m
  use twostream_bvproblem_m
  use twostream_intensity_m
  use twostream_thermalsup_m
  use twostream_l_inputs_m
  use twostream_l_writemodules_m
  use twostream_lp_miscsetups_m
  use twostream_la_solutions_m
  use twostream_lp_bvproblem_m
  use twostream_ls_bvproblem_m
  use twostream_lp_jacobians_m
  use twostream_ls_jacobians_m
  use twostream_thermalsup_plus_m
  use twostream_lssl_jacobians_m

  ! Arguments
  integer(c_int), intent(in) :: maxlayers
  integer(c_int), intent(in) :: maxtotal
  integer(c_int), intent(in) :: maxmessages
  integer(c_int), intent(in) :: maxbeams
  integer(c_int), intent(in) :: max_geometries
  integer(c_int), intent(in) :: max_user_streams
  integer(c_int), intent(in) :: max_user_relazms
  integer(c_int), intent(in) :: max_user_obsgeoms
  integer(c_int), intent(in) :: max_atmoswfs
  integer(c_int), intent(in) :: max_surfacewfs
  integer(c_int), intent(in) :: max_sleavewfs
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(in) :: do_plane_parallel
  logical(c_bool), intent(in) :: do_2s_levelout
  logical(c_bool), intent(in) :: do_mvout_only
  logical(c_bool), intent(in) :: do_additional_mvout
  logical(c_bool), intent(in) :: do_solar_sources
  logical(c_bool), intent(in) :: do_thermal_emission
  logical(c_bool), intent(in) :: do_surface_emission
  logical(c_bool), intent(in) :: do_d2s_scaling
  logical(c_bool), intent(in) :: do_brdf_surface
  logical(c_bool), intent(in) :: do_user_obsgeoms
  logical(c_bool), intent(in) :: do_surface_leaving
  logical(c_bool), intent(in) :: do_sl_isotropic
  logical(c_bool), intent(in) :: do_pentadiag_inverse
  integer(c_int), intent(in) :: bvpindex
  real(c_double), intent(in) :: bvpscalefactor
  integer(c_int), intent(in) :: taylor_order
  real(c_double), intent(in) :: taylor_small
  real(c_double), intent(in) :: tcutoff
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: ntotal
  real(c_double), intent(in) :: stream_value
  integer(c_int), intent(in) :: n_user_obsgeoms
  real(c_double), dimension(MAX_USER_OBSGEOMS, 3), intent(in) :: user_obsgeoms
  integer(c_int), intent(inout) :: n_user_streams
  real(c_double), dimension(MAX_USER_STREAMS), intent(inout) :: user_angles
  integer(c_int), intent(inout) :: n_user_relazms
  real(c_double), dimension(MAX_USER_RELAZMS), intent(inout) :: user_relazms
  real(c_double), intent(in) :: flux_factor
  integer(c_int), intent(inout) :: nbeams
  real(c_double), dimension(MAXBEAMS), intent(inout) :: beam_szas
  real(c_double), intent(inout) :: earth_radius
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: height_grid
  real(c_double), dimension(MAXLAYERS), intent(in) :: deltau_input
  real(c_double), dimension(MAXLAYERS), intent(in) :: omega_input
  real(c_double), dimension(MAXLAYERS), intent(in) :: asymm_input
  real(c_double), dimension(MAXLAYERS), intent(in) :: d2s_scaling
  real(c_double), dimension(0:MAXLAYERS), intent(in) :: thermal_bb_input
  real(c_double), intent(in) :: lambertian_albedo
  real(c_double), dimension(0:1, MAXBEAMS), intent(in) :: brdf_f_0
  real(c_double), dimension(0:1), intent(in) :: brdf_f
  real(c_double), dimension(0:1, MAX_USER_STREAMS), intent(in) :: ubrdf_f
  real(c_double), intent(in) :: emissivity
  real(c_double), intent(in) :: surfbb
  real(c_double), dimension(MAXBEAMS), intent(in) :: slterm_isotropic
  real(c_double), dimension(0:1, MAXBEAMS), intent(in) :: slterm_f_0
  logical(c_bool), intent(in) :: do_profile_wfs
  logical(c_bool), intent(in) :: do_surface_wfs
  logical(c_bool), intent(in) :: do_sleave_wfs
  logical(c_bool), dimension(MAXLAYERS), intent(in) :: layer_vary_flag
  integer(c_int), dimension(MAXLAYERS), intent(in) :: layer_vary_number
  integer(c_int), intent(in) :: n_surface_wfs
  integer(c_int), intent(in) :: n_sleave_wfs
  real(c_double), dimension(MAX_SLEAVEWFS, MAXBEAMS), intent(in) :: lssl_slterm_isotropic
  real(c_double), dimension(MAX_SLEAVEWFS, 0:1, MAXBEAMS), intent(in) :: lssl_slterm_f_0
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_deltau_input
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_omega_input
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_asymm_input
  real(c_double), dimension(MAXLAYERS, MAX_ATMOSWFS), intent(in) :: l_d2s_scaling
  real(c_double), dimension(MAX_SURFACEWFS, 0:1, MAXBEAMS), intent(in) :: ls_brdf_f_0
  real(c_double), dimension(MAX_SURFACEWFS, 0:1), intent(in) :: ls_brdf_f
  real(c_double), dimension(MAX_SURFACEWFS, 0:1, MAX_USER_STREAMS), intent(in) :: ls_ubrdf_f
  real(c_double), dimension(MAX_SURFACEWFS), intent(in) :: ls_emissivity
  real(c_double), dimension(MAX_GEOMETRIES), intent(inout) :: intensity_toa
  real(c_double), dimension(MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profilewf_toa
  real(c_double), dimension(MAX_GEOMETRIES, MAX_SURFACEWFS), intent(inout) :: surfacewf_toa
  real(c_double), dimension(MAX_GEOMETRIES), intent(inout) :: intensity_boa
  real(c_double), dimension(MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profilewf_boa
  real(c_double), dimension(MAX_GEOMETRIES, MAX_SURFACEWFS), intent(inout) :: surfacewf_boa
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS), intent(inout) :: radlevel_up
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS), intent(inout) :: radlevel_dn
  integer(c_int), intent(inout) :: n_geometries
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profjaclevel_up
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profjaclevel_dn
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS), intent(inout) :: surfjaclevel_up
  real(c_double), dimension(MAX_GEOMETRIES, 0:MAXLAYERS, MAX_SURFACEWFS), intent(inout) :: surfjaclevel_dn
  real(c_double), dimension(MAXBEAMS, 2), intent(inout) :: fluxes_toa
  real(c_double), dimension(MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profjacfluxes_toa
  real(c_double), dimension(MAXBEAMS, 2, MAX_SURFACEWFS), intent(inout) :: surfjacfluxes_toa
  real(c_double), dimension(MAXBEAMS, 2), intent(inout) :: fluxes_boa
  real(c_double), dimension(MAXBEAMS, 2, MAXLAYERS, MAX_ATMOSWFS), intent(inout) :: profjacfluxes_boa
  real(c_double), dimension(MAXBEAMS, 2, MAX_SURFACEWFS), intent(inout) :: surfjacfluxes_boa
  integer(c_int), intent(out) :: status_inputcheck
  integer(c_int), intent(out) :: c_nmessages
  integer(c_int), intent(in) :: c_messages_shape_1
  integer(c_int), intent(in) :: c_messages_len
  character(kind=c_char) , intent(inout) :: c_messages(c_messages_shape_1, c_messages_len+1)
  integer(c_int), intent(in) :: c_actions_shape_1
  integer(c_int), intent(in) :: c_actions_len
  character(kind=c_char) , intent(inout) :: c_actions(c_actions_shape_1, c_actions_len+1)
  integer(c_int), intent(out) :: status_execution
  integer(c_int), intent(in) :: e_message_len
  character(kind=c_char) , intent(inout) :: e_message(e_message_len+1)
  integer(c_int), intent(in) :: e_trace_1_len
  character(kind=c_char) , intent(inout) :: e_trace_1(e_trace_1_len+1)
  integer(c_int), intent(in) :: e_trace_2_len
  character(kind=c_char) , intent(inout) :: e_trace_2(e_trace_2_len+1)

  ! Local variables
  logical(kind=4) :: do_upwelling_lcl
  logical(kind=4) :: do_dnwelling_lcl
  logical(kind=4) :: do_plane_parallel_lcl
  logical(kind=4) :: do_2s_levelout_lcl
  logical(kind=4) :: do_mvout_only_lcl
  logical(kind=4) :: do_additional_mvout_lcl
  logical(kind=4) :: do_solar_sources_lcl
  logical(kind=4) :: do_thermal_emission_lcl
  logical(kind=4) :: do_surface_emission_lcl
  logical(kind=4) :: do_d2s_scaling_lcl
  logical(kind=4) :: do_brdf_surface_lcl
  logical(kind=4) :: do_user_obsgeoms_lcl
  logical(kind=4) :: do_surface_leaving_lcl
  logical(kind=4) :: do_sl_isotropic_lcl
  logical(kind=4) :: do_pentadiag_inverse_lcl
  logical(kind=4) :: do_profile_wfs_lcl
  logical(kind=4) :: do_surface_wfs_lcl
  logical(kind=4) :: do_sleave_wfs_lcl
  logical(kind=4), dimension(MAXLAYERS) :: layer_vary_flag_lcl
  character(kind=c_char, len=c_messages_len), dimension(0:MAXMESSAGES) :: c_messages_lcl
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx
  character(kind=c_char, len=c_actions_len), dimension(0:MAXMESSAGES) :: c_actions_lcl
  character(kind=c_char, len=e_message_len) :: e_message_lcl
  character(kind=c_char, len=e_trace_1_len) :: e_trace_1_lcl
  character(kind=c_char, len=e_trace_2_len) :: e_trace_2_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_plane_parallel_lcl = do_plane_parallel
  do_2s_levelout_lcl = do_2s_levelout
  do_mvout_only_lcl = do_mvout_only
  do_additional_mvout_lcl = do_additional_mvout
  do_solar_sources_lcl = do_solar_sources
  do_thermal_emission_lcl = do_thermal_emission
  do_surface_emission_lcl = do_surface_emission
  do_d2s_scaling_lcl = do_d2s_scaling
  do_brdf_surface_lcl = do_brdf_surface
  do_user_obsgeoms_lcl = do_user_obsgeoms
  do_surface_leaving_lcl = do_surface_leaving
  do_sl_isotropic_lcl = do_sl_isotropic
  do_pentadiag_inverse_lcl = do_pentadiag_inverse
  do_profile_wfs_lcl = do_profile_wfs
  do_surface_wfs_lcl = do_surface_wfs
  do_sleave_wfs_lcl = do_sleave_wfs
  layer_vary_flag_lcl = layer_vary_flag

  call twostream_lps_master(maxlayers, &
                            maxtotal, &
                            maxmessages, &
                            maxbeams, &
                            max_geometries, &
                            max_user_streams, &
                            max_user_relazms, &
                            max_user_obsgeoms, &
                            max_atmoswfs, &
                            max_surfacewfs, &
                            max_sleavewfs, &
                            do_upwelling_lcl, &
                            do_dnwelling_lcl, &
                            do_plane_parallel_lcl, &
                            do_2s_levelout_lcl, &
                            do_mvout_only_lcl, &
                            do_additional_mvout_lcl, &
                            do_solar_sources_lcl, &
                            do_thermal_emission_lcl, &
                            do_surface_emission_lcl, &
                            do_d2s_scaling_lcl, &
                            do_brdf_surface_lcl, &
                            do_user_obsgeoms_lcl, &
                            do_surface_leaving_lcl, &
                            do_sl_isotropic_lcl, &
                            do_pentadiag_inverse_lcl, &
                            bvpindex, &
                            bvpscalefactor, &
                            taylor_order, &
                            taylor_small, &
                            tcutoff, &
                            nlayers, &
                            ntotal, &
                            stream_value, &
                            n_user_obsgeoms, &
                            user_obsgeoms, &
                            n_user_streams, &
                            user_angles, &
                            n_user_relazms, &
                            user_relazms, &
                            flux_factor, &
                            nbeams, &
                            beam_szas, &
                            earth_radius, &
                            height_grid, &
                            deltau_input, &
                            omega_input, &
                            asymm_input, &
                            d2s_scaling, &
                            thermal_bb_input, &
                            lambertian_albedo, &
                            brdf_f_0, &
                            brdf_f, &
                            ubrdf_f, &
                            emissivity, &
                            surfbb, &
                            slterm_isotropic, &
                            slterm_f_0, &
                            do_profile_wfs_lcl, &
                            do_surface_wfs_lcl, &
                            do_sleave_wfs_lcl, &
                            layer_vary_flag_lcl, &
                            layer_vary_number, &
                            n_surface_wfs, &
                            n_sleave_wfs, &
                            lssl_slterm_isotropic, &
                            lssl_slterm_f_0, &
                            l_deltau_input, &
                            l_omega_input, &
                            l_asymm_input, &
                            l_d2s_scaling, &
                            ls_brdf_f_0, &
                            ls_brdf_f, &
                            ls_ubrdf_f, &
                            ls_emissivity, &
                            intensity_toa, &
                            profilewf_toa, &
                            surfacewf_toa, &
                            intensity_boa, &
                            profilewf_boa, &
                            surfacewf_boa, &
                            radlevel_up, &
                            radlevel_dn, &
                            n_geometries, &
                            profjaclevel_up, &
                            profjaclevel_dn, &
                            surfjaclevel_up, &
                            surfjaclevel_dn, &
                            fluxes_toa, &
                            profjacfluxes_toa, &
                            surfjacfluxes_toa, &
                            fluxes_boa, &
                            profjacfluxes_boa, &
                            surfjacfluxes_boa, &
                            status_inputcheck, &
                            c_nmessages, &
                            c_messages_lcl, &
                            c_actions_lcl, &
                            status_execution, &
                            e_message_lcl, &
                            e_trace_1_lcl, &
                            e_trace_2_lcl)

  ! Convert output arguments
  lb_1 = lbound(c_messages_lcl,1)
  do dim_idx_1 = 1, min(c_nmessages, c_messages_shape_1)
    do len_idx = 1, c_messages_len
      c_messages(dim_idx_1, len_idx) = &
          c_messages_lcl(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(c_messages_lcl(dim_idx_1-1+lb_1)(:))+1
    c_messages(dim_idx_1, len_idx) = c_null_char
  end do
  lb_1 = lbound(c_actions_lcl,1)
  do dim_idx_1 = 1, min(c_nmessages, c_actions_shape_1)
    do len_idx = 1, c_actions_len
      c_actions(dim_idx_1, len_idx) = &
          c_actions_lcl(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(c_actions_lcl(dim_idx_1-1+lb_1)(:))+1
    c_actions(dim_idx_1, len_idx) = c_null_char
  end do
  do len_idx = 1, e_message_len
    e_message(len_idx) = &
      e_message_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(e_message_lcl(:))+1
  e_message(len_idx) = c_null_char
  do len_idx = 1, e_trace_1_len
    e_trace_1(len_idx) = &
      e_trace_1_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(e_trace_1_lcl(:))+1
  e_trace_1(len_idx) = c_null_char
  do len_idx = 1, e_trace_2_len
    e_trace_2(len_idx) = &
      e_trace_2_lcl(len_idx:len_idx)
  end do
  len_idx = len_trim(e_trace_2_lcl(:))+1
  e_trace_2(len_idx) = c_null_char

end subroutine twostream_lps_master_m_twostream_lps_master_wrap

end module TWOSTREAM_LPS_MASTER_M_WRAP
