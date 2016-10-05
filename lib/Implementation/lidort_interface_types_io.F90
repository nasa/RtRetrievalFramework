module lidort_interface_types_io

use iso_c_binding
use lidort_pars

! This module was auto-generated 

implicit none

contains

! Links to type: "brdf_linsup_inputs" from module: "brdf_linsup_inputs_def" in file: "brdf_lin_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_inputs_c_write(lun, fortran_type_c) bind(C)
  use brdf_linsup_inputs_def, only : brdf_linsup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_linsup_inputs_f_write(lun, fortran_type_f)

end subroutine brdf_linsup_inputs_c_write

subroutine brdf_linsup_inputs_f_write(lun, fortran_type_f) 
  use brdf_linsup_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_linsup_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_do_kernel_factor_wfs
  write(UNIT=lun) fortran_type_f%bs_do_kernel_params_wfs
  write(UNIT=lun) fortran_type_f%bs_do_kparams_derivs
  write(UNIT=lun) fortran_type_f%bs_n_surface_wfs
  write(UNIT=lun) fortran_type_f%bs_n_kernel_factor_wfs
  write(UNIT=lun) fortran_type_f%bs_n_kernel_params_wfs
  
end subroutine brdf_linsup_inputs_f_write

subroutine brdf_linsup_inputs_c_read(lun, fortran_type_c) bind(C)
  use brdf_linsup_inputs_def, only : brdf_linsup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_linsup_inputs_f_read(lun, fortran_type_f)

end subroutine brdf_linsup_inputs_c_read

subroutine brdf_linsup_inputs_f_read(lun, fortran_type_f) 
  use brdf_linsup_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_linsup_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_do_kernel_factor_wfs
  read(UNIT=lun) fortran_type_f%bs_do_kernel_params_wfs
  read(UNIT=lun) fortran_type_f%bs_do_kparams_derivs
  read(UNIT=lun) fortran_type_f%bs_n_surface_wfs
  read(UNIT=lun) fortran_type_f%bs_n_kernel_factor_wfs
  read(UNIT=lun) fortran_type_f%bs_n_kernel_params_wfs
  
end subroutine brdf_linsup_inputs_f_read

! Links to type: "brdf_linsup_outputs" from module: "brdf_linsup_outputs_def" in file: "brdf_lin_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_linsup_outputs_c_write(lun, fortran_type_c) bind(C)
  use brdf_linsup_outputs_def, only : brdf_linsup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_linsup_outputs_f_write(lun, fortran_type_f)

end subroutine brdf_linsup_outputs_c_write

subroutine brdf_linsup_outputs_f_write(lun, fortran_type_f) 
  use brdf_linsup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_linsup_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_ls_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%bs_ls_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_ls_brdf_f
  write(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f
  write(UNIT=lun) fortran_type_f%bs_ls_emissivity
  write(UNIT=lun) fortran_type_f%bs_ls_user_emissivity
  
end subroutine brdf_linsup_outputs_f_write

subroutine brdf_linsup_outputs_c_read(lun, fortran_type_c) bind(C)
  use brdf_linsup_outputs_def, only : brdf_linsup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_linsup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_linsup_outputs_f_read(lun, fortran_type_f)

end subroutine brdf_linsup_outputs_c_read

subroutine brdf_linsup_outputs_f_read(lun, fortran_type_f) 
  use brdf_linsup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_linsup_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_ls_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%bs_ls_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_ls_brdf_f
  read(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_ls_user_brdf_f
  read(UNIT=lun) fortran_type_f%bs_ls_emissivity
  read(UNIT=lun) fortran_type_f%bs_ls_user_emissivity
  
end subroutine brdf_linsup_outputs_f_read

! Links to type: "brdf_sup_inputs" from module: "brdf_sup_inputs_def" in file: "brdf_sup_inputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_inputs_c_write(lun, fortran_type_c) bind(C)
  use brdf_sup_inputs_def, only : brdf_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_sup_inputs_f_write(lun, fortran_type_f)

end subroutine brdf_sup_inputs_c_write

subroutine brdf_sup_inputs_f_write(lun, fortran_type_f) 
  use brdf_sup_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_sup_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_do_user_streams
  write(UNIT=lun) fortran_type_f%bs_do_brdf_surface
  write(UNIT=lun) fortran_type_f%bs_do_surface_emission
  write(UNIT=lun) fortran_type_f%bs_nstreams
  write(UNIT=lun) fortran_type_f%bs_nbeams
  write(UNIT=lun) fortran_type_f%bs_beam_szas
  write(UNIT=lun) fortran_type_f%bs_n_user_relazms
  write(UNIT=lun) fortran_type_f%bs_user_relazms
  write(UNIT=lun) fortran_type_f%bs_n_user_streams
  write(UNIT=lun) fortran_type_f%bs_user_angles_input
  write(UNIT=lun) fortran_type_f%bs_n_brdf_kernels
  write(UNIT=lun) fortran_type_f%bs_brdf_names
  write(UNIT=lun) fortran_type_f%bs_which_brdf
  write(UNIT=lun) fortran_type_f%bs_n_brdf_parameters
  write(UNIT=lun) fortran_type_f%bs_brdf_parameters
  write(UNIT=lun) fortran_type_f%bs_lambertian_kernel_flag
  write(UNIT=lun) fortran_type_f%bs_brdf_factors
  write(UNIT=lun) fortran_type_f%bs_nstreams_brdf
  write(UNIT=lun) fortran_type_f%bs_do_shadow_effect
  write(UNIT=lun) fortran_type_f%bs_do_exactonly
  write(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr
  write(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr_exactonly
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_order
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nmuquad
  write(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nphiquad
  
end subroutine brdf_sup_inputs_f_write

subroutine brdf_sup_inputs_c_read(lun, fortran_type_c) bind(C)
  use brdf_sup_inputs_def, only : brdf_sup_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_sup_inputs_f_read(lun, fortran_type_f)

end subroutine brdf_sup_inputs_c_read

subroutine brdf_sup_inputs_f_read(lun, fortran_type_f) 
  use brdf_sup_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_sup_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_do_user_streams
  read(UNIT=lun) fortran_type_f%bs_do_brdf_surface
  read(UNIT=lun) fortran_type_f%bs_do_surface_emission
  read(UNIT=lun) fortran_type_f%bs_nstreams
  read(UNIT=lun) fortran_type_f%bs_nbeams
  read(UNIT=lun) fortran_type_f%bs_beam_szas
  read(UNIT=lun) fortran_type_f%bs_n_user_relazms
  read(UNIT=lun) fortran_type_f%bs_user_relazms
  read(UNIT=lun) fortran_type_f%bs_n_user_streams
  read(UNIT=lun) fortran_type_f%bs_user_angles_input
  read(UNIT=lun) fortran_type_f%bs_n_brdf_kernels
  read(UNIT=lun) fortran_type_f%bs_brdf_names
  read(UNIT=lun) fortran_type_f%bs_which_brdf
  read(UNIT=lun) fortran_type_f%bs_n_brdf_parameters
  read(UNIT=lun) fortran_type_f%bs_brdf_parameters
  read(UNIT=lun) fortran_type_f%bs_lambertian_kernel_flag
  read(UNIT=lun) fortran_type_f%bs_brdf_factors
  read(UNIT=lun) fortran_type_f%bs_nstreams_brdf
  read(UNIT=lun) fortran_type_f%bs_do_shadow_effect
  read(UNIT=lun) fortran_type_f%bs_do_exactonly
  read(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr
  read(UNIT=lun) fortran_type_f%bs_do_glitter_msrcorr_exactonly
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_order
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nmuquad
  read(UNIT=lun) fortran_type_f%bs_glitter_msrcorr_nphiquad
  
end subroutine brdf_sup_inputs_f_read

! Links to type: "brdf_sup_outputs" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_sup_outputs_c_write(lun, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_sup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_sup_outputs_f_write(lun, fortran_type_f)

end subroutine brdf_sup_outputs_c_write

subroutine brdf_sup_outputs_f_write(lun, fortran_type_f) 
  use brdf_sup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_sup_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%bs_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_brdf_f
  write(UNIT=lun) fortran_type_f%bs_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%bs_user_brdf_f
  write(UNIT=lun) fortran_type_f%bs_emissivity
  write(UNIT=lun) fortran_type_f%bs_user_emissivity
  
end subroutine brdf_sup_outputs_f_write

subroutine brdf_sup_outputs_c_read(lun, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_sup_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_sup_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_sup_outputs_f_read(lun, fortran_type_f)

end subroutine brdf_sup_outputs_c_read

subroutine brdf_sup_outputs_f_read(lun, fortran_type_f) 
  use brdf_sup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_sup_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%bs_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_brdf_f
  read(UNIT=lun) fortran_type_f%bs_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%bs_user_brdf_f
  read(UNIT=lun) fortran_type_f%bs_emissivity
  read(UNIT=lun) fortran_type_f%bs_user_emissivity
  
end subroutine brdf_sup_outputs_f_read

! Links to type: "brdf_input_exception_handling" from module: "brdf_sup_outputs_def" in file: "brdf_sup_outputs_def.F90"
! Allocs and initializes type
subroutine brdf_input_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_input_exception_handling_f_write(lun, fortran_type_f)

end subroutine brdf_input_exception_handling_c_write

subroutine brdf_input_exception_handling_f_write(lun, fortran_type_f) 
  use brdf_sup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_input_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%bs_status_inputread
  write(UNIT=lun) fortran_type_f%bs_ninputmessages
  write(UNIT=lun) fortran_type_f%bs_inputmessages
  write(UNIT=lun) fortran_type_f%bs_inputactions
  
end subroutine brdf_input_exception_handling_f_write

subroutine brdf_input_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use brdf_sup_outputs_def, only : brdf_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(brdf_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call brdf_input_exception_handling_f_read(lun, fortran_type_f)

end subroutine brdf_input_exception_handling_c_read

subroutine brdf_input_exception_handling_f_read(lun, fortran_type_f) 
  use brdf_sup_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(brdf_input_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%bs_status_inputread
  read(UNIT=lun) fortran_type_f%bs_ninputmessages
  read(UNIT=lun) fortran_type_f%bs_inputmessages
  read(UNIT=lun) fortran_type_f%bs_inputactions
  
end subroutine brdf_input_exception_handling_f_read

! Links to type: "lidort_fixed_lincontrol" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lincontrol_c_write(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_lincontrol_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_lincontrol_c_write

subroutine lidort_fixed_lincontrol_f_write(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_lincontrol), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_column_linearization
  write(UNIT=lun) fortran_type_f%ts_do_profile_linearization
  write(UNIT=lun) fortran_type_f%ts_do_surface_linearization
  write(UNIT=lun) fortran_type_f%ts_do_sleave_wfs
  write(UNIT=lun) fortran_type_f%ts_layer_vary_flag
  write(UNIT=lun) fortran_type_f%ts_layer_vary_number
  write(UNIT=lun) fortran_type_f%ts_n_totalcolumn_wfs
  write(UNIT=lun) fortran_type_f%ts_n_surface_wfs
  write(UNIT=lun) fortran_type_f%ts_n_sleave_wfs
  
end subroutine lidort_fixed_lincontrol_f_write

subroutine lidort_fixed_lincontrol_c_read(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lincontrol

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lincontrol), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_lincontrol_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_lincontrol_c_read

subroutine lidort_fixed_lincontrol_f_read(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_lincontrol), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_column_linearization
  read(UNIT=lun) fortran_type_f%ts_do_profile_linearization
  read(UNIT=lun) fortran_type_f%ts_do_surface_linearization
  read(UNIT=lun) fortran_type_f%ts_do_sleave_wfs
  read(UNIT=lun) fortran_type_f%ts_layer_vary_flag
  read(UNIT=lun) fortran_type_f%ts_layer_vary_number
  read(UNIT=lun) fortran_type_f%ts_n_totalcolumn_wfs
  read(UNIT=lun) fortran_type_f%ts_n_surface_wfs
  read(UNIT=lun) fortran_type_f%ts_n_sleave_wfs
  
end subroutine lidort_fixed_lincontrol_f_read

! Links to type: "lidort_fixed_linoptical" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_linoptical_c_write(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_linoptical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_linoptical_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_linoptical_c_write

subroutine lidort_fixed_linoptical_f_write(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_linoptical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_l_deltau_vert_input
  write(UNIT=lun) fortran_type_f%ts_l_omega_total_input
  write(UNIT=lun) fortran_type_f%ts_l_phasmoms_total_input
  
end subroutine lidort_fixed_linoptical_f_write

subroutine lidort_fixed_linoptical_c_read(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_linoptical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_linoptical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_linoptical_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_linoptical_c_read

subroutine lidort_fixed_linoptical_f_read(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_linoptical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_l_deltau_vert_input
  read(UNIT=lun) fortran_type_f%ts_l_omega_total_input
  read(UNIT=lun) fortran_type_f%ts_l_phasmoms_total_input
  
end subroutine lidort_fixed_linoptical_f_read

! Links to type: "lidort_fixed_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_lininputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_lininputs_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_lininputs_c_write

subroutine lidort_fixed_lininputs_f_write(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_lininputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_fixed_lincontrol), pointer :: cont_lcl  
  type(lidort_fixed_linoptical), pointer :: optical_lcl  
  
  ! Get pointer to types
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  
  call lidort_fixed_lincontrol_f_write(lun, cont_lcl)
  call lidort_fixed_linoptical_f_write(lun, optical_lcl)
  
end subroutine lidort_fixed_lininputs_f_write

subroutine lidort_fixed_lininputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_fixed_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_lininputs_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_lininputs_c_read

subroutine lidort_fixed_lininputs_f_read(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_lininputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_fixed_lincontrol), pointer :: cont_lcl  
  type(lidort_fixed_linoptical), pointer :: optical_lcl  
  
  ! Get pointer to types
  cont_lcl => fortran_type_f%cont
  optical_lcl => fortran_type_f%optical
  
  call lidort_fixed_lincontrol_f_read(lun, cont_lcl)
  call lidort_fixed_linoptical_f_read(lun, optical_lcl)
  
end subroutine lidort_fixed_lininputs_f_read

! Links to type: "lidort_modified_lininputs" from module: "lidort_lininputs_def" in file: "lidort_lin_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_lininputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_modified_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_lininputs_f_write(lun, fortran_type_f)

end subroutine lidort_modified_lininputs_c_write

subroutine lidort_modified_lininputs_f_write(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_lininputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%dummy
  
end subroutine lidort_modified_lininputs_f_write

subroutine lidort_modified_lininputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_lininputs_def, only : lidort_modified_lininputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_lininputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_lininputs_f_read(lun, fortran_type_f)

end subroutine lidort_modified_lininputs_c_read

subroutine lidort_modified_lininputs_f_read(lun, fortran_type_f) 
  use lidort_lininputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_lininputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%dummy
  
end subroutine lidort_modified_lininputs_f_read

! Links to type: "lidort_linatmos" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linatmos_c_write(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linatmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linatmos_f_write(lun, fortran_type_f)

end subroutine lidort_linatmos_c_write

subroutine lidort_linatmos_f_write(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linatmos), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_columnwf
  write(UNIT=lun) fortran_type_f%ts_mint_columnwf
  write(UNIT=lun) fortran_type_f%ts_flux_columnwf
  write(UNIT=lun) fortran_type_f%ts_profilewf
  write(UNIT=lun) fortran_type_f%ts_mint_profilewf
  write(UNIT=lun) fortran_type_f%ts_flux_profilewf
  
end subroutine lidort_linatmos_f_write

subroutine lidort_linatmos_c_read(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linatmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linatmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linatmos_f_read(lun, fortran_type_f)

end subroutine lidort_linatmos_c_read

subroutine lidort_linatmos_f_read(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linatmos), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_columnwf
  read(UNIT=lun) fortran_type_f%ts_mint_columnwf
  read(UNIT=lun) fortran_type_f%ts_flux_columnwf
  read(UNIT=lun) fortran_type_f%ts_profilewf
  read(UNIT=lun) fortran_type_f%ts_mint_profilewf
  read(UNIT=lun) fortran_type_f%ts_flux_profilewf
  
end subroutine lidort_linatmos_f_read

! Links to type: "lidort_linsurf" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linsurf_c_write(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linsurf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsurf_f_write(lun, fortran_type_f)

end subroutine lidort_linsurf_c_write

subroutine lidort_linsurf_f_write(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsurf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_surfacewf
  write(UNIT=lun) fortran_type_f%ts_mint_surfacewf
  write(UNIT=lun) fortran_type_f%ts_flux_surfacewf
  
end subroutine lidort_linsurf_f_write

subroutine lidort_linsurf_c_read(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linsurf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsurf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsurf_f_read(lun, fortran_type_f)

end subroutine lidort_linsurf_c_read

subroutine lidort_linsurf_f_read(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsurf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_surfacewf
  read(UNIT=lun) fortran_type_f%ts_mint_surfacewf
  read(UNIT=lun) fortran_type_f%ts_flux_surfacewf
  
end subroutine lidort_linsurf_f_read

! Links to type: "lidort_linoutputs" from module: "lidort_linoutputs_def" in file: "lidort_lin_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_linoutputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linoutputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linoutputs_f_write(lun, fortran_type_f)

end subroutine lidort_linoutputs_c_write

subroutine lidort_linoutputs_f_write(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linoutputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linatmos), pointer :: atmos_lcl  
  type(lidort_linsurf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  surf_lcl => fortran_type_f%surf
  
  call lidort_linatmos_f_write(lun, atmos_lcl)
  call lidort_linsurf_f_write(lun, surf_lcl)
  
end subroutine lidort_linoutputs_f_write

subroutine lidort_linoutputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_linoutputs_def, only : lidort_linoutputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linoutputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linoutputs_f_read(lun, fortran_type_f)

end subroutine lidort_linoutputs_c_read

subroutine lidort_linoutputs_f_read(lun, fortran_type_f) 
  use lidort_linoutputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linoutputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linatmos), pointer :: atmos_lcl  
  type(lidort_linsurf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  surf_lcl => fortran_type_f%surf
  
  call lidort_linatmos_f_read(lun, atmos_lcl)
  call lidort_linsurf_f_read(lun, surf_lcl)
  
end subroutine lidort_linoutputs_f_read

! Links to type: "lidort_linsup_brdf" from module: "lidort_linsup_brdf_def" in file: "lidort_lin_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_brdf_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_brdf_def, only : lidort_linsup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_brdf_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_brdf_c_write

subroutine lidort_linsup_brdf_f_write(lun, fortran_type_f) 
  use lidort_linsup_brdf_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_brdf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_ls_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%ts_ls_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_ls_brdf_f
  write(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f
  write(UNIT=lun) fortran_type_f%ts_ls_emissivity
  write(UNIT=lun) fortran_type_f%ts_ls_user_emissivity
  
end subroutine lidort_linsup_brdf_f_write

subroutine lidort_linsup_brdf_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_brdf_def, only : lidort_linsup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_brdf_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_brdf_c_read

subroutine lidort_linsup_brdf_f_read(lun, fortran_type_f) 
  use lidort_linsup_brdf_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_brdf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_ls_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%ts_ls_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_ls_brdf_f
  read(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_ls_user_brdf_f
  read(UNIT=lun) fortran_type_f%ts_ls_emissivity
  read(UNIT=lun) fortran_type_f%ts_ls_user_emissivity
  
end subroutine lidort_linsup_brdf_f_read

! Links to type: "lidort_linsup_ss_atmos" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_atmos_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_atmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_atmos_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_ss_atmos_c_write

subroutine lidort_linsup_ss_atmos_f_write(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss_atmos), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_columnwf_ss
  write(UNIT=lun) fortran_type_f%ts_columnwf_db
  write(UNIT=lun) fortran_type_f%ts_profilewf_ss
  write(UNIT=lun) fortran_type_f%ts_profilewf_db
  
end subroutine lidort_linsup_ss_atmos_f_write

subroutine lidort_linsup_ss_atmos_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_atmos

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_atmos), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_atmos_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_ss_atmos_c_read

subroutine lidort_linsup_ss_atmos_f_read(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss_atmos), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_columnwf_ss
  read(UNIT=lun) fortran_type_f%ts_columnwf_db
  read(UNIT=lun) fortran_type_f%ts_profilewf_ss
  read(UNIT=lun) fortran_type_f%ts_profilewf_db
  
end subroutine lidort_linsup_ss_atmos_f_read

! Links to type: "lidort_linsup_ss_surf" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_surf_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_surf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_surf_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_ss_surf_c_write

subroutine lidort_linsup_ss_surf_f_write(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss_surf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_surfacewf_db
  
end subroutine lidort_linsup_ss_surf_f_write

subroutine lidort_linsup_ss_surf_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss_surf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss_surf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_surf_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_ss_surf_c_read

subroutine lidort_linsup_ss_surf_f_read(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss_surf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_surfacewf_db
  
end subroutine lidort_linsup_ss_surf_f_read

! Links to type: "lidort_linsup_ss" from module: "lidort_linsup_ss_def" in file: "lidort_lin_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_ss_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_ss_c_write

subroutine lidort_linsup_ss_f_write(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linsup_ss_atmos), pointer :: atmos_lcl  
  type(lidort_linsup_ss_surf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  surf_lcl => fortran_type_f%surf
  
  call lidort_linsup_ss_atmos_f_write(lun, atmos_lcl)
  call lidort_linsup_ss_surf_f_write(lun, surf_lcl)
  
end subroutine lidort_linsup_ss_f_write

subroutine lidort_linsup_ss_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_ss_def, only : lidort_linsup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_ss_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_ss_c_read

subroutine lidort_linsup_ss_f_read(lun, fortran_type_f) 
  use lidort_linsup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_ss), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linsup_ss_atmos), pointer :: atmos_lcl  
  type(lidort_linsup_ss_surf), pointer :: surf_lcl  
  
  ! Get pointer to types
  atmos_lcl => fortran_type_f%atmos
  surf_lcl => fortran_type_f%surf
  
  call lidort_linsup_ss_atmos_f_read(lun, atmos_lcl)
  call lidort_linsup_ss_surf_f_read(lun, surf_lcl)
  
end subroutine lidort_linsup_ss_f_read

! Links to type: "lidort_linsup_sleave" from module: "lidort_linsup_sleave_def" in file: "lidort_lin_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_sleave_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_sleave_def, only : lidort_linsup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_sleave_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_sleave_c_write

subroutine lidort_linsup_sleave_f_write(lun, fortran_type_f) 
  use lidort_linsup_sleave_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_sleave), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_isotropic
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_userangles
  write(UNIT=lun) fortran_type_f%ts_lssl_slterm_f_0
  write(UNIT=lun) fortran_type_f%ts_lssl_user_slterm_f_0
  
end subroutine lidort_linsup_sleave_f_write

subroutine lidort_linsup_sleave_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_sleave_def, only : lidort_linsup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_sleave_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_sleave_c_read

subroutine lidort_linsup_sleave_f_read(lun, fortran_type_f) 
  use lidort_linsup_sleave_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_linsup_sleave), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_isotropic
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_userangles
  read(UNIT=lun) fortran_type_f%ts_lssl_slterm_f_0
  read(UNIT=lun) fortran_type_f%ts_lssl_user_slterm_f_0
  
end subroutine lidort_linsup_sleave_f_read

! Links to type: "lidort_linsup_inout" from module: "lidort_linsup_inout_def" in file: "lidort_lin_sup_def.F90"
! Allocs and initializes type
subroutine lidort_linsup_inout_c_write(lun, fortran_type_c) bind(C)
  use lidort_linsup_inout_def, only : lidort_linsup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_inout_f_write(lun, fortran_type_f)

end subroutine lidort_linsup_inout_c_write

subroutine lidort_linsup_inout_f_write(lun, fortran_type_f) 
  use lidort_linsup_inout_def
  use lidort_linsup_brdf_def
  use lidort_linsup_ss_def
  use lidort_linsup_sleave_def
  
  integer, intent(in) :: lun
  type(lidort_linsup_inout), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linsup_brdf), pointer :: brdf_lcl  
  type(lidort_linsup_ss), pointer :: ss_lcl  
  type(lidort_linsup_sleave), pointer :: sleave_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  ss_lcl => fortran_type_f%ss
  sleave_lcl => fortran_type_f%sleave
  
  call lidort_linsup_brdf_f_write(lun, brdf_lcl)
  call lidort_linsup_ss_f_write(lun, ss_lcl)
  call lidort_linsup_sleave_f_write(lun, sleave_lcl)
  
end subroutine lidort_linsup_inout_f_write

subroutine lidort_linsup_inout_c_read(lun, fortran_type_c) bind(C)
  use lidort_linsup_inout_def, only : lidort_linsup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_linsup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_linsup_inout_f_read(lun, fortran_type_f)

end subroutine lidort_linsup_inout_c_read

subroutine lidort_linsup_inout_f_read(lun, fortran_type_f) 
  use lidort_linsup_inout_def
  use lidort_linsup_brdf_def
  use lidort_linsup_ss_def
  use lidort_linsup_sleave_def
  
  integer, intent(in) :: lun
  type(lidort_linsup_inout), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_linsup_brdf), pointer :: brdf_lcl  
  type(lidort_linsup_ss), pointer :: ss_lcl  
  type(lidort_linsup_sleave), pointer :: sleave_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  ss_lcl => fortran_type_f%ss
  sleave_lcl => fortran_type_f%sleave
  
  call lidort_linsup_brdf_f_read(lun, brdf_lcl)
  call lidort_linsup_ss_f_read(lun, ss_lcl)
  call lidort_linsup_sleave_f_read(lun, sleave_lcl)
  
end subroutine lidort_linsup_inout_f_read

! Links to type: "lidort_main_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_main_outputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_main_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_main_outputs_f_write(lun, fortran_type_f)

end subroutine lidort_main_outputs_c_write

subroutine lidort_main_outputs_f_write(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_main_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_intensity
  write(UNIT=lun) fortran_type_f%ts_mean_intensity
  write(UNIT=lun) fortran_type_f%ts_flux_integral
  write(UNIT=lun) fortran_type_f%ts_dnflux_direct
  write(UNIT=lun) fortran_type_f%ts_dnmean_direct
  write(UNIT=lun) fortran_type_f%ts_fourier_saved
  write(UNIT=lun) fortran_type_f%ts_n_geometries
  
end subroutine lidort_main_outputs_f_write

subroutine lidort_main_outputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_main_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_main_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_main_outputs_f_read(lun, fortran_type_f)

end subroutine lidort_main_outputs_c_read

subroutine lidort_main_outputs_f_read(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_main_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_intensity
  read(UNIT=lun) fortran_type_f%ts_mean_intensity
  read(UNIT=lun) fortran_type_f%ts_flux_integral
  read(UNIT=lun) fortran_type_f%ts_dnflux_direct
  read(UNIT=lun) fortran_type_f%ts_dnmean_direct
  read(UNIT=lun) fortran_type_f%ts_fourier_saved
  read(UNIT=lun) fortran_type_f%ts_n_geometries
  
end subroutine lidort_main_outputs_f_read

! Links to type: "lidort_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_exception_handling_f_write(lun, fortran_type_f)

end subroutine lidort_exception_handling_c_write

subroutine lidort_exception_handling_f_write(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_status_inputcheck
  write(UNIT=lun) fortran_type_f%ts_ncheckmessages
  write(UNIT=lun) fortran_type_f%ts_checkmessages
  write(UNIT=lun) fortran_type_f%ts_actions
  write(UNIT=lun) fortran_type_f%ts_status_calculation
  write(UNIT=lun) fortran_type_f%ts_message
  write(UNIT=lun) fortran_type_f%ts_trace_1
  write(UNIT=lun) fortran_type_f%ts_trace_2
  write(UNIT=lun) fortran_type_f%ts_trace_3
  
end subroutine lidort_exception_handling_f_write

subroutine lidort_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_exception_handling_f_read(lun, fortran_type_f)

end subroutine lidort_exception_handling_c_read

subroutine lidort_exception_handling_f_read(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_status_inputcheck
  read(UNIT=lun) fortran_type_f%ts_ncheckmessages
  read(UNIT=lun) fortran_type_f%ts_checkmessages
  read(UNIT=lun) fortran_type_f%ts_actions
  read(UNIT=lun) fortran_type_f%ts_status_calculation
  read(UNIT=lun) fortran_type_f%ts_message
  read(UNIT=lun) fortran_type_f%ts_trace_1
  read(UNIT=lun) fortran_type_f%ts_trace_2
  read(UNIT=lun) fortran_type_f%ts_trace_3
  
end subroutine lidort_exception_handling_f_read

! Links to type: "lidort_input_exception_handling" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_input_exception_handling_c_write(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_input_exception_handling_f_write(lun, fortran_type_f)

end subroutine lidort_input_exception_handling_c_write

subroutine lidort_input_exception_handling_f_write(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_input_exception_handling), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_status_inputread
  write(UNIT=lun) fortran_type_f%ts_ninputmessages
  write(UNIT=lun) fortran_type_f%ts_inputmessages
  write(UNIT=lun) fortran_type_f%ts_inputactions
  
end subroutine lidort_input_exception_handling_f_write

subroutine lidort_input_exception_handling_c_read(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_input_exception_handling

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_input_exception_handling), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_input_exception_handling_f_read(lun, fortran_type_f)

end subroutine lidort_input_exception_handling_c_read

subroutine lidort_input_exception_handling_f_read(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_input_exception_handling), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_status_inputread
  read(UNIT=lun) fortran_type_f%ts_ninputmessages
  read(UNIT=lun) fortran_type_f%ts_inputmessages
  read(UNIT=lun) fortran_type_f%ts_inputactions
  
end subroutine lidort_input_exception_handling_f_read

! Links to type: "lidort_outputs" from module: "lidort_outputs_def" in file: "lidort_outputs_def.F90"
! Allocs and initializes type
subroutine lidort_outputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_outputs_f_write(lun, fortran_type_f)

end subroutine lidort_outputs_c_write

subroutine lidort_outputs_f_write(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_outputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_main_outputs), pointer :: main_lcl  
  type(lidort_exception_handling), pointer :: status_lcl  
  
  ! Get pointer to types
  main_lcl => fortran_type_f%main
  status_lcl => fortran_type_f%status
  
  call lidort_main_outputs_f_write(lun, main_lcl)
  call lidort_exception_handling_f_write(lun, status_lcl)
  
end subroutine lidort_outputs_f_write

subroutine lidort_outputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_outputs_def, only : lidort_outputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_outputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_outputs_f_read(lun, fortran_type_f)

end subroutine lidort_outputs_c_read

subroutine lidort_outputs_f_read(lun, fortran_type_f) 
  use lidort_outputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_outputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_main_outputs), pointer :: main_lcl  
  type(lidort_exception_handling), pointer :: status_lcl  
  
  ! Get pointer to types
  main_lcl => fortran_type_f%main
  status_lcl => fortran_type_f%status
  
  call lidort_main_outputs_f_read(lun, main_lcl)
  call lidort_exception_handling_f_read(lun, status_lcl)
  
end subroutine lidort_outputs_f_read

! Links to type: "lidort_sup_brdf" from module: "lidort_sup_brdf_def" in file: "lidort_sup_brdf_def.F90"
! Allocs and initializes type
subroutine lidort_sup_brdf_c_write(lun, fortran_type_c) bind(C)
  use lidort_sup_brdf_def, only : lidort_sup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_brdf_f_write(lun, fortran_type_f)

end subroutine lidort_sup_brdf_c_write

subroutine lidort_sup_brdf_f_write(lun, fortran_type_f) 
  use lidort_sup_brdf_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_brdf), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_exactdb_brdfunc
  write(UNIT=lun) fortran_type_f%ts_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_brdf_f
  write(UNIT=lun) fortran_type_f%ts_user_brdf_f_0
  write(UNIT=lun) fortran_type_f%ts_user_brdf_f
  write(UNIT=lun) fortran_type_f%ts_emissivity
  write(UNIT=lun) fortran_type_f%ts_user_emissivity
  
end subroutine lidort_sup_brdf_f_write

subroutine lidort_sup_brdf_c_read(lun, fortran_type_c) bind(C)
  use lidort_sup_brdf_def, only : lidort_sup_brdf

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_brdf), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_brdf_f_read(lun, fortran_type_f)

end subroutine lidort_sup_brdf_c_read

subroutine lidort_sup_brdf_f_read(lun, fortran_type_f) 
  use lidort_sup_brdf_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_brdf), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_exactdb_brdfunc
  read(UNIT=lun) fortran_type_f%ts_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_brdf_f
  read(UNIT=lun) fortran_type_f%ts_user_brdf_f_0
  read(UNIT=lun) fortran_type_f%ts_user_brdf_f
  read(UNIT=lun) fortran_type_f%ts_emissivity
  read(UNIT=lun) fortran_type_f%ts_user_emissivity
  
end subroutine lidort_sup_brdf_f_read

! Links to type: "lidort_sup_sleave" from module: "lidort_sup_sleave_def" in file: "lidort_sup_sleave_def.F90"
! Allocs and initializes type
subroutine lidort_sup_sleave_c_write(lun, fortran_type_c) bind(C)
  use lidort_sup_sleave_def, only : lidort_sup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_sleave_f_write(lun, fortran_type_f)

end subroutine lidort_sup_sleave_c_write

subroutine lidort_sup_sleave_f_write(lun, fortran_type_f) 
  use lidort_sup_sleave_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_sleave), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_slterm_isotropic
  write(UNIT=lun) fortran_type_f%ts_slterm_userangles
  write(UNIT=lun) fortran_type_f%ts_slterm_f_0
  write(UNIT=lun) fortran_type_f%ts_user_slterm_f_0
  
end subroutine lidort_sup_sleave_f_write

subroutine lidort_sup_sleave_c_read(lun, fortran_type_c) bind(C)
  use lidort_sup_sleave_def, only : lidort_sup_sleave

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_sleave), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_sleave_f_read(lun, fortran_type_f)

end subroutine lidort_sup_sleave_c_read

subroutine lidort_sup_sleave_f_read(lun, fortran_type_f) 
  use lidort_sup_sleave_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_sleave), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_slterm_isotropic
  read(UNIT=lun) fortran_type_f%ts_slterm_userangles
  read(UNIT=lun) fortran_type_f%ts_slterm_f_0
  read(UNIT=lun) fortran_type_f%ts_user_slterm_f_0
  
end subroutine lidort_sup_sleave_f_read

! Links to type: "lidort_sup_ss" from module: "lidort_sup_ss_def" in file: "lidort_sup_ss_def.F90"
! Allocs and initializes type
subroutine lidort_sup_ss_c_write(lun, fortran_type_c) bind(C)
  use lidort_sup_ss_def, only : lidort_sup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_ss_f_write(lun, fortran_type_f)

end subroutine lidort_sup_ss_c_write

subroutine lidort_sup_ss_f_write(lun, fortran_type_f) 
  use lidort_sup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_ss), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_intensity_ss
  write(UNIT=lun) fortran_type_f%ts_intensity_db
  
end subroutine lidort_sup_ss_f_write

subroutine lidort_sup_ss_c_read(lun, fortran_type_c) bind(C)
  use lidort_sup_ss_def, only : lidort_sup_ss

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_ss), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_ss_f_read(lun, fortran_type_f)

end subroutine lidort_sup_ss_c_read

subroutine lidort_sup_ss_f_read(lun, fortran_type_f) 
  use lidort_sup_ss_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_sup_ss), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_intensity_ss
  read(UNIT=lun) fortran_type_f%ts_intensity_db
  
end subroutine lidort_sup_ss_f_read

! Links to type: "lidort_sup_inout" from module: "lidort_sup_inout_def" in file: "lidort_sup_def.F90"
! Allocs and initializes type
subroutine lidort_sup_inout_c_write(lun, fortran_type_c) bind(C)
  use lidort_sup_inout_def, only : lidort_sup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_inout_f_write(lun, fortran_type_f)

end subroutine lidort_sup_inout_c_write

subroutine lidort_sup_inout_f_write(lun, fortran_type_f) 
  use lidort_sup_inout_def
  use lidort_sup_brdf_def
  use lidort_sup_ss_def
  use lidort_sup_sleave_def
  
  integer, intent(in) :: lun
  type(lidort_sup_inout), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_sup_brdf), pointer :: brdf_lcl  
  type(lidort_sup_ss), pointer :: ss_lcl  
  type(lidort_sup_sleave), pointer :: sleave_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  ss_lcl => fortran_type_f%ss
  sleave_lcl => fortran_type_f%sleave
  
  call lidort_sup_brdf_f_write(lun, brdf_lcl)
  call lidort_sup_ss_f_write(lun, ss_lcl)
  call lidort_sup_sleave_f_write(lun, sleave_lcl)
  
end subroutine lidort_sup_inout_f_write

subroutine lidort_sup_inout_c_read(lun, fortran_type_c) bind(C)
  use lidort_sup_inout_def, only : lidort_sup_inout

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_sup_inout), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_sup_inout_f_read(lun, fortran_type_f)

end subroutine lidort_sup_inout_c_read

subroutine lidort_sup_inout_f_read(lun, fortran_type_f) 
  use lidort_sup_inout_def
  use lidort_sup_brdf_def
  use lidort_sup_ss_def
  use lidort_sup_sleave_def
  
  integer, intent(in) :: lun
  type(lidort_sup_inout), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_sup_brdf), pointer :: brdf_lcl  
  type(lidort_sup_ss), pointer :: ss_lcl  
  type(lidort_sup_sleave), pointer :: sleave_lcl  
  
  ! Get pointer to types
  brdf_lcl => fortran_type_f%brdf
  ss_lcl => fortran_type_f%ss
  sleave_lcl => fortran_type_f%sleave
  
  call lidort_sup_brdf_f_read(lun, brdf_lcl)
  call lidort_sup_ss_f_read(lun, ss_lcl)
  call lidort_sup_sleave_f_read(lun, sleave_lcl)
  
end subroutine lidort_sup_inout_f_read

! Links to type: "lidort_fixed_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_boolean_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_boolean_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_boolean_c_write

subroutine lidort_fixed_boolean_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_boolean), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_fullrad_mode
  write(UNIT=lun) fortran_type_f%ts_do_sscorr_truncation
  write(UNIT=lun) fortran_type_f%ts_do_ss_external
  write(UNIT=lun) fortran_type_f%ts_do_ssfull
  write(UNIT=lun) fortran_type_f%ts_do_thermal_emission
  write(UNIT=lun) fortran_type_f%ts_do_surface_emission
  write(UNIT=lun) fortran_type_f%ts_do_plane_parallel
  write(UNIT=lun) fortran_type_f%ts_do_brdf_surface
  write(UNIT=lun) fortran_type_f%ts_do_upwelling
  write(UNIT=lun) fortran_type_f%ts_do_dnwelling
  write(UNIT=lun) fortran_type_f%ts_do_surface_leaving
  write(UNIT=lun) fortran_type_f%ts_do_sl_isotropic
  
end subroutine lidort_fixed_boolean_f_write

subroutine lidort_fixed_boolean_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_boolean_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_boolean_c_read

subroutine lidort_fixed_boolean_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_boolean), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_fullrad_mode
  read(UNIT=lun) fortran_type_f%ts_do_sscorr_truncation
  read(UNIT=lun) fortran_type_f%ts_do_ss_external
  read(UNIT=lun) fortran_type_f%ts_do_ssfull
  read(UNIT=lun) fortran_type_f%ts_do_thermal_emission
  read(UNIT=lun) fortran_type_f%ts_do_surface_emission
  read(UNIT=lun) fortran_type_f%ts_do_plane_parallel
  read(UNIT=lun) fortran_type_f%ts_do_brdf_surface
  read(UNIT=lun) fortran_type_f%ts_do_upwelling
  read(UNIT=lun) fortran_type_f%ts_do_dnwelling
  read(UNIT=lun) fortran_type_f%ts_do_surface_leaving
  read(UNIT=lun) fortran_type_f%ts_do_sl_isotropic
  
end subroutine lidort_fixed_boolean_f_read

! Links to type: "lidort_fixed_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_control_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_control_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_control_c_write

subroutine lidort_fixed_control_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_control), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_nstreams
  write(UNIT=lun) fortran_type_f%ts_nlayers
  write(UNIT=lun) fortran_type_f%ts_nfinelayers
  write(UNIT=lun) fortran_type_f%ts_n_thermal_coeffs
  write(UNIT=lun) fortran_type_f%ts_lidort_accuracy
  
end subroutine lidort_fixed_control_f_write

subroutine lidort_fixed_control_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_control_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_control_c_read

subroutine lidort_fixed_control_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_control), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_nstreams
  read(UNIT=lun) fortran_type_f%ts_nlayers
  read(UNIT=lun) fortran_type_f%ts_nfinelayers
  read(UNIT=lun) fortran_type_f%ts_n_thermal_coeffs
  read(UNIT=lun) fortran_type_f%ts_lidort_accuracy
  
end subroutine lidort_fixed_control_f_read

! Links to type: "lidort_fixed_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_sunrays_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_sunrays_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_sunrays_c_write

subroutine lidort_fixed_sunrays_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_sunrays), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_flux_factor
  
end subroutine lidort_fixed_sunrays_f_write

subroutine lidort_fixed_sunrays_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_sunrays_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_sunrays_c_read

subroutine lidort_fixed_sunrays_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_sunrays), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_flux_factor
  
end subroutine lidort_fixed_sunrays_f_read

! Links to type: "lidort_fixed_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_uservalues_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_uservalues_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_uservalues_c_write

subroutine lidort_fixed_uservalues_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_uservalues), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_n_user_streams
  write(UNIT=lun) fortran_type_f%ts_n_user_levels
  
end subroutine lidort_fixed_uservalues_f_write

subroutine lidort_fixed_uservalues_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_uservalues_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_uservalues_c_read

subroutine lidort_fixed_uservalues_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_uservalues), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_n_user_streams
  read(UNIT=lun) fortran_type_f%ts_n_user_levels
  
end subroutine lidort_fixed_uservalues_f_read

! Links to type: "lidort_fixed_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_chapman_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_chapman_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_chapman_c_write

subroutine lidort_fixed_chapman_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_chapman), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_height_grid
  write(UNIT=lun) fortran_type_f%ts_pressure_grid
  write(UNIT=lun) fortran_type_f%ts_temperature_grid
  write(UNIT=lun) fortran_type_f%ts_finegrid
  write(UNIT=lun) fortran_type_f%ts_rfindex_parameter
  
end subroutine lidort_fixed_chapman_f_write

subroutine lidort_fixed_chapman_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_chapman_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_chapman_c_read

subroutine lidort_fixed_chapman_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_chapman), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_height_grid
  read(UNIT=lun) fortran_type_f%ts_pressure_grid
  read(UNIT=lun) fortran_type_f%ts_temperature_grid
  read(UNIT=lun) fortran_type_f%ts_finegrid
  read(UNIT=lun) fortran_type_f%ts_rfindex_parameter
  
end subroutine lidort_fixed_chapman_f_read

! Links to type: "lidort_fixed_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_optical_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_optical_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_optical_c_write

subroutine lidort_fixed_optical_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_optical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_deltau_vert_input
  write(UNIT=lun) fortran_type_f%ts_phasmoms_total_input
  write(UNIT=lun) fortran_type_f%ts_thermal_bb_input
  write(UNIT=lun) fortran_type_f%ts_lambertian_albedo
  write(UNIT=lun) fortran_type_f%ts_surface_bb_input
  
end subroutine lidort_fixed_optical_f_write

subroutine lidort_fixed_optical_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_optical_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_optical_c_read

subroutine lidort_fixed_optical_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_optical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_deltau_vert_input
  read(UNIT=lun) fortran_type_f%ts_phasmoms_total_input
  read(UNIT=lun) fortran_type_f%ts_thermal_bb_input
  read(UNIT=lun) fortran_type_f%ts_lambertian_albedo
  read(UNIT=lun) fortran_type_f%ts_surface_bb_input
  
end subroutine lidort_fixed_optical_f_read

! Links to type: "lidort_fixed_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_fixed_inputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_inputs_f_write(lun, fortran_type_f)

end subroutine lidort_fixed_inputs_c_write

subroutine lidort_fixed_inputs_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_fixed_boolean), pointer :: bool_lcl  
  type(lidort_fixed_control), pointer :: cont_lcl  
  type(lidort_fixed_sunrays), pointer :: sunrays_lcl  
  type(lidort_fixed_uservalues), pointer :: userval_lcl  
  type(lidort_fixed_chapman), pointer :: chapman_lcl  
  type(lidort_fixed_optical), pointer :: optical_lcl  
  
  ! Get pointer to types
  bool_lcl => fortran_type_f%bool
  cont_lcl => fortran_type_f%cont
  sunrays_lcl => fortran_type_f%sunrays
  userval_lcl => fortran_type_f%userval
  chapman_lcl => fortran_type_f%chapman
  optical_lcl => fortran_type_f%optical
  
  call lidort_fixed_boolean_f_write(lun, bool_lcl)
  call lidort_fixed_control_f_write(lun, cont_lcl)
  call lidort_fixed_sunrays_f_write(lun, sunrays_lcl)
  call lidort_fixed_uservalues_f_write(lun, userval_lcl)
  call lidort_fixed_chapman_f_write(lun, chapman_lcl)
  call lidort_fixed_optical_f_write(lun, optical_lcl)
  
end subroutine lidort_fixed_inputs_f_write

subroutine lidort_fixed_inputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_fixed_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_fixed_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_fixed_inputs_f_read(lun, fortran_type_f)

end subroutine lidort_fixed_inputs_c_read

subroutine lidort_fixed_inputs_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_fixed_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_fixed_boolean), pointer :: bool_lcl  
  type(lidort_fixed_control), pointer :: cont_lcl  
  type(lidort_fixed_sunrays), pointer :: sunrays_lcl  
  type(lidort_fixed_uservalues), pointer :: userval_lcl  
  type(lidort_fixed_chapman), pointer :: chapman_lcl  
  type(lidort_fixed_optical), pointer :: optical_lcl  
  
  ! Get pointer to types
  bool_lcl => fortran_type_f%bool
  cont_lcl => fortran_type_f%cont
  sunrays_lcl => fortran_type_f%sunrays
  userval_lcl => fortran_type_f%userval
  chapman_lcl => fortran_type_f%chapman
  optical_lcl => fortran_type_f%optical
  
  call lidort_fixed_boolean_f_read(lun, bool_lcl)
  call lidort_fixed_control_f_read(lun, cont_lcl)
  call lidort_fixed_sunrays_f_read(lun, sunrays_lcl)
  call lidort_fixed_uservalues_f_read(lun, userval_lcl)
  call lidort_fixed_chapman_f_read(lun, chapman_lcl)
  call lidort_fixed_optical_f_read(lun, optical_lcl)
  
end subroutine lidort_fixed_inputs_f_read

! Links to type: "lidort_modified_boolean" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_boolean_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_boolean_f_write(lun, fortran_type_f)

end subroutine lidort_modified_boolean_c_write

subroutine lidort_modified_boolean_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_boolean), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_do_sscorr_nadir
  write(UNIT=lun) fortran_type_f%ts_do_sscorr_outgoing
  write(UNIT=lun) fortran_type_f%ts_do_double_convtest
  write(UNIT=lun) fortran_type_f%ts_do_solar_sources
  write(UNIT=lun) fortran_type_f%ts_do_refractive_geometry
  write(UNIT=lun) fortran_type_f%ts_do_chapman_function
  write(UNIT=lun) fortran_type_f%ts_do_rayleigh_only
  write(UNIT=lun) fortran_type_f%ts_do_isotropic_only
  write(UNIT=lun) fortran_type_f%ts_do_no_azimuth
  write(UNIT=lun) fortran_type_f%ts_do_all_fourier
  write(UNIT=lun) fortran_type_f%ts_do_deltam_scaling
  write(UNIT=lun) fortran_type_f%ts_do_solution_saving
  write(UNIT=lun) fortran_type_f%ts_do_bvp_telescoping
  write(UNIT=lun) fortran_type_f%ts_do_user_streams
  write(UNIT=lun) fortran_type_f%ts_do_additional_mvout
  write(UNIT=lun) fortran_type_f%ts_do_mvout_only
  write(UNIT=lun) fortran_type_f%ts_do_thermal_transonly
  
end subroutine lidort_modified_boolean_f_write

subroutine lidort_modified_boolean_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_boolean

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_boolean), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_boolean_f_read(lun, fortran_type_f)

end subroutine lidort_modified_boolean_c_read

subroutine lidort_modified_boolean_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_boolean), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_do_sscorr_nadir
  read(UNIT=lun) fortran_type_f%ts_do_sscorr_outgoing
  read(UNIT=lun) fortran_type_f%ts_do_double_convtest
  read(UNIT=lun) fortran_type_f%ts_do_solar_sources
  read(UNIT=lun) fortran_type_f%ts_do_refractive_geometry
  read(UNIT=lun) fortran_type_f%ts_do_chapman_function
  read(UNIT=lun) fortran_type_f%ts_do_rayleigh_only
  read(UNIT=lun) fortran_type_f%ts_do_isotropic_only
  read(UNIT=lun) fortran_type_f%ts_do_no_azimuth
  read(UNIT=lun) fortran_type_f%ts_do_all_fourier
  read(UNIT=lun) fortran_type_f%ts_do_deltam_scaling
  read(UNIT=lun) fortran_type_f%ts_do_solution_saving
  read(UNIT=lun) fortran_type_f%ts_do_bvp_telescoping
  read(UNIT=lun) fortran_type_f%ts_do_user_streams
  read(UNIT=lun) fortran_type_f%ts_do_additional_mvout
  read(UNIT=lun) fortran_type_f%ts_do_mvout_only
  read(UNIT=lun) fortran_type_f%ts_do_thermal_transonly
  
end subroutine lidort_modified_boolean_f_read

! Links to type: "lidort_modified_control" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_control_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_control_f_write(lun, fortran_type_f)

end subroutine lidort_modified_control_c_write

subroutine lidort_modified_control_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_control), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_nmoments_input
  
end subroutine lidort_modified_control_f_write

subroutine lidort_modified_control_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_control

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_control), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_control_f_read(lun, fortran_type_f)

end subroutine lidort_modified_control_c_read

subroutine lidort_modified_control_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_control), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_nmoments_input
  
end subroutine lidort_modified_control_f_read

! Links to type: "lidort_modified_sunrays" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_sunrays_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_sunrays_f_write(lun, fortran_type_f)

end subroutine lidort_modified_sunrays_c_write

subroutine lidort_modified_sunrays_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_sunrays), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_nbeams
  write(UNIT=lun) fortran_type_f%ts_beam_szas
  
end subroutine lidort_modified_sunrays_f_write

subroutine lidort_modified_sunrays_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_sunrays

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_sunrays), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_sunrays_f_read(lun, fortran_type_f)

end subroutine lidort_modified_sunrays_c_read

subroutine lidort_modified_sunrays_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_sunrays), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_nbeams
  read(UNIT=lun) fortran_type_f%ts_beam_szas
  
end subroutine lidort_modified_sunrays_f_read

! Links to type: "lidort_modified_uservalues" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_uservalues_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_uservalues_f_write(lun, fortran_type_f)

end subroutine lidort_modified_uservalues_c_write

subroutine lidort_modified_uservalues_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_uservalues), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_n_user_relazms
  write(UNIT=lun) fortran_type_f%ts_user_relazms
  write(UNIT=lun) fortran_type_f%ts_user_angles_input
  write(UNIT=lun) fortran_type_f%ts_user_levels
  write(UNIT=lun) fortran_type_f%ts_geometry_specheight
  
end subroutine lidort_modified_uservalues_f_write

subroutine lidort_modified_uservalues_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_uservalues

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_uservalues), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_uservalues_f_read(lun, fortran_type_f)

end subroutine lidort_modified_uservalues_c_read

subroutine lidort_modified_uservalues_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_uservalues), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_n_user_relazms
  read(UNIT=lun) fortran_type_f%ts_user_relazms
  read(UNIT=lun) fortran_type_f%ts_user_angles_input
  read(UNIT=lun) fortran_type_f%ts_user_levels
  read(UNIT=lun) fortran_type_f%ts_geometry_specheight
  
end subroutine lidort_modified_uservalues_f_read

! Links to type: "lidort_modified_chapman" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_chapman_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_chapman_f_write(lun, fortran_type_f)

end subroutine lidort_modified_chapman_c_write

subroutine lidort_modified_chapman_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_chapman), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_earth_radius
  
end subroutine lidort_modified_chapman_f_write

subroutine lidort_modified_chapman_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_chapman

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_chapman), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_chapman_f_read(lun, fortran_type_f)

end subroutine lidort_modified_chapman_c_read

subroutine lidort_modified_chapman_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_chapman), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_earth_radius
  
end subroutine lidort_modified_chapman_f_read

! Links to type: "lidort_modified_optical" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_optical_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_optical_f_write(lun, fortran_type_f)

end subroutine lidort_modified_optical_c_write

subroutine lidort_modified_optical_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_optical), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  write(UNIT=lun) fortran_type_f%ts_omega_total_input
  
end subroutine lidort_modified_optical_f_write

subroutine lidort_modified_optical_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_optical

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_optical), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_optical_f_read(lun, fortran_type_f)

end subroutine lidort_modified_optical_c_read

subroutine lidort_modified_optical_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_optical), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  
  ! Get pointer to types
  
  read(UNIT=lun) fortran_type_f%ts_omega_total_input
  
end subroutine lidort_modified_optical_f_read

! Links to type: "lidort_modified_inputs" from module: "lidort_inputs_def" in file: "lidort_inputs_def.F90"
! Allocs and initializes type
subroutine lidort_modified_inputs_c_write(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(in)    :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_inputs_f_write(lun, fortran_type_f)

end subroutine lidort_modified_inputs_c_write

subroutine lidort_modified_inputs_f_write(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_inputs), intent(in), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_modified_boolean), pointer :: mbool_lcl  
  type(lidort_modified_control), pointer :: mcont_lcl  
  type(lidort_modified_sunrays), pointer :: msunrays_lcl  
  type(lidort_modified_uservalues), pointer :: muserval_lcl  
  type(lidort_modified_chapman), pointer :: mchapman_lcl  
  type(lidort_modified_optical), pointer :: moptical_lcl  
  
  ! Get pointer to types
  mbool_lcl => fortran_type_f%mbool
  mcont_lcl => fortran_type_f%mcont
  msunrays_lcl => fortran_type_f%msunrays
  muserval_lcl => fortran_type_f%muserval
  mchapman_lcl => fortran_type_f%mchapman
  moptical_lcl => fortran_type_f%moptical
  
  call lidort_modified_boolean_f_write(lun, mbool_lcl)
  call lidort_modified_control_f_write(lun, mcont_lcl)
  call lidort_modified_sunrays_f_write(lun, msunrays_lcl)
  call lidort_modified_uservalues_f_write(lun, muserval_lcl)
  call lidort_modified_chapman_f_write(lun, mchapman_lcl)
  call lidort_modified_optical_f_write(lun, moptical_lcl)
  
end subroutine lidort_modified_inputs_f_write

subroutine lidort_modified_inputs_c_read(lun, fortran_type_c) bind(C)
  use lidort_inputs_def, only : lidort_modified_inputs

  integer(c_int), intent(in) :: lun
  type(c_ptr), intent(inout) :: fortran_type_c

  type(lidort_modified_inputs), pointer :: fortran_type_f

  call c_f_pointer(fortran_type_c, fortran_type_f)

  call lidort_modified_inputs_f_read(lun, fortran_type_f)

end subroutine lidort_modified_inputs_c_read

subroutine lidort_modified_inputs_f_read(lun, fortran_type_f) 
  use lidort_inputs_def
  use lidort_pars
  
  integer, intent(in) :: lun
  type(lidort_modified_inputs), intent(inout), pointer :: fortran_type_f

  ! Type pointers declarations
  type(lidort_modified_boolean), pointer :: mbool_lcl  
  type(lidort_modified_control), pointer :: mcont_lcl  
  type(lidort_modified_sunrays), pointer :: msunrays_lcl  
  type(lidort_modified_uservalues), pointer :: muserval_lcl  
  type(lidort_modified_chapman), pointer :: mchapman_lcl  
  type(lidort_modified_optical), pointer :: moptical_lcl  
  
  ! Get pointer to types
  mbool_lcl => fortran_type_f%mbool
  mcont_lcl => fortran_type_f%mcont
  msunrays_lcl => fortran_type_f%msunrays
  muserval_lcl => fortran_type_f%muserval
  mchapman_lcl => fortran_type_f%mchapman
  moptical_lcl => fortran_type_f%moptical
  
  call lidort_modified_boolean_f_read(lun, mbool_lcl)
  call lidort_modified_control_f_read(lun, mcont_lcl)
  call lidort_modified_sunrays_f_read(lun, msunrays_lcl)
  call lidort_modified_uservalues_f_read(lun, muserval_lcl)
  call lidort_modified_chapman_f_read(lun, mchapman_lcl)
  call lidort_modified_optical_f_read(lun, moptical_lcl)
  
end subroutine lidort_modified_inputs_f_read


end module lidort_interface_types_io