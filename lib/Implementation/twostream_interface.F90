
module TWOSTREAM_LS_BRDF_SUPPLEMENT_M_WRAP

use iso_c_binding
use twostream_ls_brdf_supplement_m

! This module was auto-generated 

implicit none

contains

subroutine twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap (lambertian_kernel_flag, &
                                                                        do_shadow_effect, &
                                                                        do_surface_emission, &
                                                                        nbeams, &
                                                                        n_user_streams, &
                                                                        n_user_relazms, &
                                                                        beam_szas, &
                                                                        user_angles, &
                                                                        user_relazms, &
                                                                        stream_value, &
                                                                        nstreams_brdf, &
                                                                        n_brdf_kernels, &
                                                                        which_brdf, &
                                                                        brdf_factors, &
                                                                        n_brdf_parameters, &
                                                                        brdf_parameters, &
                                                                        nspars, &
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
                                                                        ls_emissivity) bind(C)

  ! Arguments
  logical(c_bool), dimension(3), intent(in) :: lambertian_kernel_flag
  logical(c_bool), intent(in) :: do_shadow_effect
  logical(c_bool), intent(in) :: do_surface_emission
  integer(c_int), intent(in) :: nbeams
  integer(c_int), intent(in) :: n_user_streams
  integer(c_int), intent(in) :: n_user_relazms
  real(c_double), dimension(NBEAMS), intent(in) :: beam_szas
  real(c_double), dimension(N_USER_STREAMS), intent(in) :: user_angles
  real(c_double), dimension(N_USER_RELAZMS), intent(in) :: user_relazms
  real(c_double), intent(in) :: stream_value
  integer(c_int), intent(in) :: nstreams_brdf
  integer(c_int), intent(in) :: n_brdf_kernels
  integer(c_int), dimension(3), intent(in) :: which_brdf
  real(c_double), dimension(3), intent(in) :: brdf_factors
  integer(c_int), dimension(3), intent(in) :: n_brdf_parameters
  real(c_double), dimension(3, 3), intent(in) :: brdf_parameters
  integer(c_int), intent(in) :: nspars
  logical(c_bool), dimension(3), intent(in) :: do_kernel_factor_wfs
  logical(c_bool), dimension(3, 3), intent(in) :: do_kernel_params_wfs
  logical(c_bool), dimension(3), intent(out) :: do_kparams_derivs
  integer(c_int), intent(out) :: n_surface_wfs
  integer(c_int), intent(out) :: n_kernel_factor_wfs
  integer(c_int), intent(out) :: n_kernel_params_wfs
  real(c_double), dimension(0:1, NBEAMS), intent(out) :: brdf_f_0
  real(c_double), dimension(0:1), intent(out) :: brdf_f
  real(c_double), dimension(0:1, N_USER_STREAMS), intent(out) :: ubrdf_f
  real(c_double), intent(out) :: emissivity
  real(c_double), dimension(NSPARS, 0:1, NBEAMS), intent(out) :: ls_brdf_f_0
  real(c_double), dimension(NSPARS, 0:1), intent(out) :: ls_brdf_f
  real(c_double), dimension(NSPARS, 0:1, N_USER_STREAMS), intent(out) :: ls_ubrdf_f
  real(c_double), dimension(NSPARS), intent(out) :: ls_emissivity

  ! Local variables
  logical(kind=4), dimension(3) :: lambertian_kernel_flag_lcl
  logical(kind=4) :: do_shadow_effect_lcl
  logical(kind=4) :: do_surface_emission_lcl
  logical(kind=4), dimension(3) :: do_kernel_factor_wfs_lcl
  logical(kind=4), dimension(3, 3) :: do_kernel_params_wfs_lcl
  logical(kind=4), dimension(3) :: do_kparams_derivs_lcl

  ! Convert input arguments
  lambertian_kernel_flag_lcl = lambertian_kernel_flag
  do_shadow_effect_lcl = do_shadow_effect
  do_surface_emission_lcl = do_surface_emission
  do_kernel_factor_wfs_lcl = do_kernel_factor_wfs
  do_kernel_params_wfs_lcl = do_kernel_params_wfs

  call twostream_ls_brdfmaster(lambertian_kernel_flag_lcl, &
                               do_shadow_effect_lcl, &
                               do_surface_emission_lcl, &
                               nbeams, &
                               n_user_streams, &
                               n_user_relazms, &
                               beam_szas, &
                               user_angles, &
                               user_relazms, &
                               stream_value, &
                               nstreams_brdf, &
                               n_brdf_kernels, &
                               which_brdf, &
                               brdf_factors, &
                               n_brdf_parameters, &
                               brdf_parameters, &
                               nspars, &
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
                               ls_emissivity)

  ! Convert output arguments
  do_kparams_derivs = do_kparams_derivs_lcl

end subroutine twostream_ls_brdf_supplement_m_twostream_ls_brdfmaster_wrap

end module TWOSTREAM_LS_BRDF_SUPPLEMENT_M_WRAP

module TWOSTREAM_L_MASTER_M_WRAP

use iso_c_binding
use twostream_l_master_m

! This module was auto-generated 

implicit none

contains

subroutine twostream_l_master_m_twostream_l_master_wrap (thread, &
                                                         do_upwelling, &
                                                         do_dnwelling, &
                                                         do_plane_parallel, &
                                                         do_solar_sources, &
                                                         do_thermal_emission, &
                                                         do_surface_emission, &
                                                         do_d2s_scaling, &
                                                         do_brdf_surface, &
                                                         pure_nadir, &
                                                         nthreads, &
                                                         nlayers, &
                                                         ntotal, &
                                                         n_geometries, &
                                                         stream_value, &
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
                                                         do_sim_only, &
                                                         do_profile_wfs, &
                                                         do_column_wfs, &
                                                         do_surface_wfs, &
                                                         npars, &
                                                         nspars, &
                                                         layer_vary_flag, &
                                                         layer_vary_number, &
                                                         n_column_wfs, &
                                                         n_surface_wfs, &
                                                         l_deltau_input, &
                                                         l_omega_input, &
                                                         l_asymm_input, &
                                                         l_d2s_scaling, &
                                                         ls_brdf_f_0, &
                                                         ls_brdf_f, &
                                                         ls_ubrdf_f, &
                                                         ls_emiss, &
                                                         intensity_toa, &
                                                         profilewf_toa, &
                                                         columnwf_toa, &
                                                         surfacewf_toa, &
                                                         intensity_boa, &
                                                         profilewf_boa, &
                                                         columnwf_boa, &
                                                         surfacewf_boa, &
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

  ! Arguments
  integer(c_int), intent(in) :: thread
  logical(c_bool), intent(in) :: do_upwelling
  logical(c_bool), intent(in) :: do_dnwelling
  logical(c_bool), intent(in) :: do_plane_parallel
  logical(c_bool), intent(in) :: do_solar_sources
  logical(c_bool), intent(in) :: do_thermal_emission
  logical(c_bool), intent(in) :: do_surface_emission
  logical(c_bool), intent(in) :: do_d2s_scaling
  logical(c_bool), intent(in) :: do_brdf_surface
  logical(c_bool), intent(in) :: pure_nadir
  integer(c_int), intent(in) :: nthreads
  integer(c_int), intent(in) :: nlayers
  integer(c_int), intent(in) :: ntotal
  integer(c_int), intent(in) :: n_geometries
  real(c_double), intent(in) :: stream_value
  integer(c_int), intent(in) :: n_user_streams
  real(c_double), dimension(N_USER_STREAMS), intent(in) :: user_angles
  integer(c_int), intent(in) :: n_user_relazms
  real(c_double), dimension(N_USER_RELAZMS), intent(in) :: user_relazms
  real(c_double), intent(in) :: flux_factor
  integer(c_int), intent(in) :: nbeams
  real(c_double), dimension(NBEAMS), intent(in) :: beam_szas
  real(c_double), intent(inout) :: earth_radius
  real(c_double), dimension(0:NLAYERS), intent(in) :: height_grid
  real(c_double), dimension(NLAYERS, NTHREADS), intent(in) :: deltau_input
  real(c_double), dimension(NLAYERS, NTHREADS), intent(in) :: omega_input
  real(c_double), dimension(NLAYERS, NTHREADS), intent(in) :: asymm_input
  real(c_double), dimension(NLAYERS, NTHREADS), intent(in) :: d2s_scaling
  real(c_double), dimension(0:NLAYERS), intent(in) :: thermal_bb_input
  real(c_double), dimension(NTHREADS), intent(in) :: lambertian_albedo
  real(c_double), dimension(0:1, NBEAMS), intent(in) :: brdf_f_0
  real(c_double), dimension(0:1), intent(in) :: brdf_f
  real(c_double), dimension(0:1, N_USER_STREAMS), intent(in) :: ubrdf_f
  real(c_double), intent(in) :: emissivity
  real(c_double), intent(in) :: surfbb
  logical(c_bool), intent(in) :: do_sim_only
  logical(c_bool), intent(in) :: do_profile_wfs
  logical(c_bool), intent(in) :: do_column_wfs
  logical(c_bool), intent(in) :: do_surface_wfs
  integer(c_int), intent(in) :: npars
  integer(c_int), intent(in) :: nspars
  logical(c_bool), dimension(NLAYERS), intent(in) :: layer_vary_flag
  integer(c_int), dimension(NLAYERS), intent(in) :: layer_vary_number
  integer(c_int), intent(in) :: n_column_wfs
  integer(c_int), intent(in) :: n_surface_wfs
  real(c_double), dimension(NLAYERS, NPARS, NTHREADS), intent(in) :: l_deltau_input
  real(c_double), dimension(NLAYERS, NPARS, NTHREADS), intent(in) :: l_omega_input
  real(c_double), dimension(NLAYERS, NPARS, NTHREADS), intent(in) :: l_asymm_input
  real(c_double), dimension(NLAYERS, NPARS, NTHREADS), intent(in) :: l_d2s_scaling
  real(c_double), dimension(NSPARS, 0:1, NBEAMS), intent(in) :: ls_brdf_f_0
  real(c_double), dimension(NSPARS, 0:1), intent(in) :: ls_brdf_f
  real(c_double), dimension(NSPARS, 0:1, N_USER_STREAMS), intent(in) :: ls_ubrdf_f
  real(c_double), dimension(NSPARS), intent(in) :: ls_emiss
  real(c_double), dimension(N_GEOMETRIES, NTHREADS), intent(inout) :: intensity_toa
  real(c_double), dimension(N_GEOMETRIES, NLAYERS, NPARS, NTHREADS), intent(inout) :: profilewf_toa
  real(c_double), dimension(N_GEOMETRIES, NPARS, NTHREADS), intent(inout) :: columnwf_toa
  real(c_double), dimension(N_GEOMETRIES, NSPARS, NTHREADS), intent(inout) :: surfacewf_toa
  real(c_double), dimension(N_GEOMETRIES, NTHREADS), intent(inout) :: intensity_boa
  real(c_double), dimension(N_GEOMETRIES, NLAYERS, NPARS, NTHREADS), intent(inout) :: profilewf_boa
  real(c_double), dimension(N_GEOMETRIES, NPARS, NTHREADS), intent(inout) :: columnwf_boa
  real(c_double), dimension(N_GEOMETRIES, NSPARS, NTHREADS), intent(inout) :: surfacewf_boa
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
  logical(kind=4) :: do_solar_sources_lcl
  logical(kind=4) :: do_thermal_emission_lcl
  logical(kind=4) :: do_surface_emission_lcl
  logical(kind=4) :: do_d2s_scaling_lcl
  logical(kind=4) :: do_brdf_surface_lcl
  logical(kind=4) :: pure_nadir_lcl
  logical(kind=4) :: do_sim_only_lcl
  logical(kind=4) :: do_profile_wfs_lcl
  logical(kind=4) :: do_column_wfs_lcl
  logical(kind=4) :: do_surface_wfs_lcl
  logical(kind=4), dimension(NLAYERS) :: layer_vary_flag_lcl
  character(kind=c_char, len=100), dimension(100) :: c_messages_lcl
  integer :: dim_idx_1
  integer :: lb_1
  integer :: len_idx
  character(kind=c_char, len=100), dimension(100) :: c_actions_lcl
  character(kind=c_char, len=e_message_len) :: e_message_lcl
  character(kind=c_char, len=e_trace_1_len) :: e_trace_1_lcl
  character(kind=c_char, len=e_trace_2_len) :: e_trace_2_lcl

  ! Convert input arguments
  do_upwelling_lcl = do_upwelling
  do_dnwelling_lcl = do_dnwelling
  do_plane_parallel_lcl = do_plane_parallel
  do_solar_sources_lcl = do_solar_sources
  do_thermal_emission_lcl = do_thermal_emission
  do_surface_emission_lcl = do_surface_emission
  do_d2s_scaling_lcl = do_d2s_scaling
  do_brdf_surface_lcl = do_brdf_surface
  pure_nadir_lcl = pure_nadir
  do_sim_only_lcl = do_sim_only
  do_profile_wfs_lcl = do_profile_wfs
  do_column_wfs_lcl = do_column_wfs
  do_surface_wfs_lcl = do_surface_wfs
  layer_vary_flag_lcl = layer_vary_flag

  call twostream_l_master(thread, &
                          do_upwelling_lcl, &
                          do_dnwelling_lcl, &
                          do_plane_parallel_lcl, &
                          do_solar_sources_lcl, &
                          do_thermal_emission_lcl, &
                          do_surface_emission_lcl, &
                          do_d2s_scaling_lcl, &
                          do_brdf_surface_lcl, &
                          pure_nadir_lcl, &
                          nthreads, &
                          nlayers, &
                          ntotal, &
                          n_geometries, &
                          stream_value, &
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
                          do_sim_only_lcl, &
                          do_profile_wfs_lcl, &
                          do_column_wfs_lcl, &
                          do_surface_wfs_lcl, &
                          npars, &
                          nspars, &
                          layer_vary_flag_lcl, &
                          layer_vary_number, &
                          n_column_wfs, &
                          n_surface_wfs, &
                          l_deltau_input, &
                          l_omega_input, &
                          l_asymm_input, &
                          l_d2s_scaling, &
                          ls_brdf_f_0, &
                          ls_brdf_f, &
                          ls_ubrdf_f, &
                          ls_emiss, &
                          intensity_toa, &
                          profilewf_toa, &
                          columnwf_toa, &
                          surfacewf_toa, &
                          intensity_boa, &
                          profilewf_boa, &
                          columnwf_boa, &
                          surfacewf_boa, &
                          status_inputcheck, &
                          c_nmessages, &
                          c_messages_lcl, &
                          c_actions_lcl, &
                          status_execution, &
                          e_message_lcl, &
                          e_trace_1_lcl, &
                          e_trace_2_lcl)

  ! Convert output arguments
  lb_1 = lbound(c_messages,1)
  do dim_idx_1 = 1, min(c_nmessages, c_messages_shape_1)
    do len_idx = 1, c_messages_len
      c_messages(dim_idx_1, len_idx) = &
          c_messages_lcl(dim_idx_1-1+lb_1)(len_idx:len_idx)
    end do
    len_idx = len_trim(c_messages_lcl(dim_idx_1-1+lb_1)(:))+1
    c_messages(dim_idx_1, len_idx) = c_null_char
  end do
  lb_1 = lbound(c_actions,1)
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

end subroutine twostream_l_master_m_twostream_l_master_wrap

end module TWOSTREAM_L_MASTER_M_WRAP