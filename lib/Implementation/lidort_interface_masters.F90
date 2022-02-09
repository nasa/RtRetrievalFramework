
module BRDF_LIN_SUP_MASTERS_M_WRAP

use iso_c_binding
use brdf_lin_sup_masters_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine brdf_lin_sup_masters_m_brdf_lin_inputmaster_wrap (filnam_in_len, &
                                                             filnam_in, &
                                                             brdf_sup_in_in, &
                                                             brdf_linsup_in_in, &
                                                             brdf_sup_inputstatus_in) bind(C)
  use lidort_pars_m
  use brdf_findpar_m
  use brdf_sup_inputs_def_m
  use brdf_sup_outputs_def_m
  use brdf_lin_sup_inputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_linsup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_linsup_inputs), pointer :: brdf_linsup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_linsup_in_in, brdf_linsup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)

  call brdf_lin_inputmaster(filnam_lcl, &
                            brdf_sup_in_lcl, &
                            brdf_linsup_in_lcl, &
                            brdf_sup_inputstatus_lcl)

end subroutine brdf_lin_sup_masters_m_brdf_lin_inputmaster_wrap

subroutine brdf_lin_sup_masters_m_brdf_lin_mainmaster_wrap (do_debug_restoration_in, &
                                                            nmoments_input_in, &
                                                            brdf_sup_in_in, &
                                                            brdf_linsup_in_in, &
                                                            brdf_sup_out_in, &
                                                            brdf_linsup_out_in, &
                                                            brdf_sup_outputstatus_in) bind(C)
  use lidort_pars_m
  use brdf_sup_inputs_def_m
  use brdf_sup_outputs_def_m
  use brdf_lin_sup_inputs_def_m
  use brdf_lin_sup_outputs_def_m
  use brdf_sup_aux_m
  use brdf_sup_kernels_m
  use brdf_lin_sup_kernels_m
  use brdf_sup_routines_m
  use brdf_lin_sup_routines_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(in) :: brdf_linsup_in_in
  type(c_ptr), intent(out) :: brdf_sup_out_in
  type(c_ptr), intent(out) :: brdf_linsup_out_in
  type(c_ptr), intent(out) :: brdf_sup_outputstatus_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_linsup_inputs), pointer :: brdf_linsup_in_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(brdf_linsup_outputs), pointer :: brdf_linsup_out_lcl
  type(brdf_output_exception_handling), pointer :: brdf_sup_outputstatus_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_linsup_in_in, brdf_linsup_in_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(brdf_linsup_out_in, brdf_linsup_out_lcl)
  call c_f_pointer(brdf_sup_outputstatus_in, brdf_sup_outputstatus_lcl)

  call brdf_lin_mainmaster(do_debug_restoration_lcl, &
                           nmoments_input_in, &
                           brdf_sup_in_lcl, &
                           brdf_linsup_in_lcl, &
                           brdf_sup_out_lcl, &
                           brdf_linsup_out_lcl, &
                           brdf_sup_outputstatus_lcl)

end subroutine brdf_lin_sup_masters_m_brdf_lin_mainmaster_wrap

end module BRDF_LIN_SUP_MASTERS_M_WRAP

module BRDF_SUP_MASTERS_M_WRAP

use iso_c_binding
use brdf_sup_masters_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine brdf_sup_masters_m_brdf_inputmaster_wrap (filnam_in_len, &
                                                     filnam_in, &
                                                     brdf_sup_in_in, &
                                                     brdf_sup_inputstatus_in) bind(C)
  use lidort_pars_m
  use brdf_findpar_m
  use brdf_sup_inputs_def_m
  use brdf_sup_outputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)

  call brdf_inputmaster(filnam_lcl, &
                        brdf_sup_in_lcl, &
                        brdf_sup_inputstatus_lcl)

end subroutine brdf_sup_masters_m_brdf_inputmaster_wrap

subroutine brdf_sup_masters_m_brdf_mainmaster_wrap (do_debug_restoration_in, &
                                                    nmoments_input_in, &
                                                    brdf_sup_in_in, &
                                                    brdf_sup_out_in, &
                                                    brdf_sup_outputstatus_in) bind(C)
  use lidort_pars_m
  use brdf_sup_inputs_def_m
  use brdf_sup_outputs_def_m
  use brdf_sup_aux_m
  use brdf_sup_kernels_m
  use brdf_sup_routines_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_sup_out_in
  type(c_ptr), intent(out) :: brdf_sup_outputstatus_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(brdf_output_exception_handling), pointer :: brdf_sup_outputstatus_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(brdf_sup_outputstatus_in, brdf_sup_outputstatus_lcl)

  call brdf_mainmaster(do_debug_restoration_lcl, &
                       nmoments_input_in, &
                       brdf_sup_in_lcl, &
                       brdf_sup_out_lcl, &
                       brdf_sup_outputstatus_lcl)

end subroutine brdf_sup_masters_m_brdf_mainmaster_wrap

end module BRDF_SUP_MASTERS_M_WRAP

module LIDORT_INPUTS_M_WRAP

use iso_c_binding
use lidort_inputs_m
use lidort_pars_m
use lidort_pars_m
use lidort_aux_m

! This module was auto-generated 

implicit none

contains

subroutine inputs_m_brdf_sup_init_wrap (lidort_sup_in) bind(C)
  use lidort_pars_m
  use lidort_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_sup_in

  ! Local variables
  type(lidort_sup_inout), pointer :: lidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)

  call lidort_brdf_sup_init(lidort_sup_lcl)

end subroutine inputs_m_brdf_sup_init_wrap

subroutine inputs_m_input_master_wrap (filnam_in_len, &
                                       filnam_in, &
                                       lidort_fixin_in, &
                                       lidort_modin_in, &
                                       lidort_inputstatus_in) bind(C)
  use lidort_pars_m
  use lidort_inputs_def_m
  use lidort_outputs_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: lidort_fixin_in
  type(c_ptr), intent(out) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_input_exception_handling), pointer :: lidort_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_inputstatus_in, lidort_inputstatus_lcl)

  call lidort_input_master(filnam_lcl, &
                           lidort_fixin_lcl, &
                           lidort_modin_lcl, &
                           lidort_inputstatus_lcl)

end subroutine inputs_m_input_master_wrap

subroutine inputs_m_sleave_sup_init_wrap (lidort_sup_in) bind(C)
  use lidort_pars_m
  use lidort_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_sup_in

  ! Local variables
  type(lidort_sup_inout), pointer :: lidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)

  call lidort_sleave_sup_init(lidort_sup_lcl)

end subroutine inputs_m_sleave_sup_init_wrap

subroutine inputs_m_ss_sup_init_wrap (lidort_sup_in) bind(C)
  use lidort_pars_m
  use lidort_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_sup_in

  ! Local variables
  type(lidort_sup_inout), pointer :: lidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)

  call lidort_ss_sup_init(lidort_sup_lcl)

end subroutine inputs_m_ss_sup_init_wrap

subroutine inputs_m_sup_init_wrap (lidort_sup_in) bind(C)
  use lidort_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_sup_in

  ! Local variables
  type(lidort_sup_inout), pointer :: lidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)

  call lidort_sup_init(lidort_sup_lcl)

end subroutine inputs_m_sup_init_wrap

end module LIDORT_INPUTS_M_WRAP

module LIDORT_MASTERS_M_WRAP

use iso_c_binding
use lidort_masters_m
use lidort_pars_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine masters_m_master_wrap (do_debug_input_in, &
                                  lidort_fixin_in, &
                                  lidort_modin_in, &
                                  lidort_sup_in, &
                                  lidort_out_in) bind(C)
  use lidort_pars_m
  use lidort_io_defs_m
  use lidort_inputs_m
  use lidort_geometry_m
  use lidort_miscsetups_m
  use lidort_bvproblem_m
  use lidort_converge_m
  use lidort_thermalsup_m
  use lidort_writemodules_m
  use lidort_sfo_interface_m
  use lidort_pars_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(out) :: lidort_out_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)

  call lidort_master(do_debug_input_lcl, &
                     lidort_fixin_lcl, &
                     lidort_modin_lcl, &
                     lidort_sup_lcl, &
                     lidort_out_lcl)

end subroutine masters_m_master_wrap

end module LIDORT_MASTERS_M_WRAP

module LIDORT_L_INPUTS_M_WRAP

use iso_c_binding
use lidort_l_inputs_m
use lidort_pars_m
use lidort_pars_m
use lidort_aux_m

! This module was auto-generated 

implicit none

contains

subroutine l_inputs_m_brdf_linsup_init_wrap (lidort_linsup_in) bind(C)
  use lidort_pars_m
  use lidort_lin_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_linsup_in

  ! Local variables
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)

  call lidort_brdf_linsup_init(lidort_linsup_lcl)

end subroutine l_inputs_m_brdf_linsup_init_wrap

subroutine l_inputs_m_l_input_master_wrap (filnam_in_len, &
                                           filnam_in, &
                                           lidort_fixin_in, &
                                           lidort_modin_in, &
                                           lidort_linfixin_in, &
                                           lidort_linmodin_in, &
                                           lidort_inputstatus_in) bind(C)
  use lidort_pars_m
  use lidort_inputs_def_m
  use lidort_lin_inputs_def_m
  use lidort_outputs_def_m
  use lidort_inputs_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  integer(c_int), intent(in) :: filnam_in_len
  character(kind=c_char) , intent(inout) :: filnam_in(filnam_in_len+1)
  type(c_ptr), intent(out) :: lidort_fixin_in
  type(c_ptr), intent(out) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_linfixin_in
  type(c_ptr), intent(out) :: lidort_linmodin_in
  type(c_ptr), intent(out) :: lidort_inputstatus_in

  ! Local variables
  character(kind=c_char, len=filnam_in_len) :: filnam_lcl
  integer :: len_idx
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_fixed_lininputs), pointer :: lidort_linfixin_lcl
  type(lidort_modified_lininputs), pointer :: lidort_linmodin_lcl
  type(lidort_input_exception_handling), pointer :: lidort_inputstatus_lcl

  ! Convert input arguments
  do len_idx = 1, filnam_in_len
    filnam_lcl(len_idx:len_idx) = &
      filnam_in(len_idx)
  end do
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_inputstatus_in, lidort_inputstatus_lcl)

  call lidort_l_input_master(filnam_lcl, &
                             lidort_fixin_lcl, &
                             lidort_modin_lcl, &
                             lidort_linfixin_lcl, &
                             lidort_linmodin_lcl, &
                             lidort_inputstatus_lcl)

end subroutine l_inputs_m_l_input_master_wrap

subroutine l_inputs_m_linsup_init_wrap (lidort_linsup_in) bind(C)
  use lidort_lin_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_linsup_in

  ! Local variables
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)

  call lidort_linsup_init(lidort_linsup_lcl)

end subroutine l_inputs_m_linsup_init_wrap

subroutine l_inputs_m_sleave_linsup_init_wrap (lidort_linsup_in) bind(C)
  use lidort_pars_m
  use lidort_lin_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_linsup_in

  ! Local variables
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)

  call lidort_sleave_linsup_init(lidort_linsup_lcl)

end subroutine l_inputs_m_sleave_linsup_init_wrap

subroutine l_inputs_m_ss_linsup_init_wrap (lidort_linsup_in) bind(C)
  use lidort_pars_m
  use lidort_lin_sup_inout_def_m
  use lidort_pars_m
  use lidort_aux_m

  ! Arguments
  type(c_ptr), intent(inout) :: lidort_linsup_in

  ! Local variables
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)

  call lidort_ss_linsup_init(lidort_linsup_lcl)

end subroutine l_inputs_m_ss_linsup_init_wrap

end module LIDORT_L_INPUTS_M_WRAP

module LIDORT_LCS_MASTERS_M_WRAP

use iso_c_binding
use lidort_lcs_masters_m
use lidort_pars_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine lcs_masters_m_lcs_master_wrap (do_debug_input_in, &
                                          lidort_fixin_in, &
                                          lidort_modin_in, &
                                          lidort_sup_in, &
                                          lidort_out_in, &
                                          lidort_linfixin_in, &
                                          lidort_linmodin_in, &
                                          lidort_linsup_in, &
                                          lidort_linout_in) bind(C)
  use lidort_pars_m
  use lidort_io_defs_m
  use lidort_lin_io_defs_m
  use lidort_inputs_m
  use lidort_l_inputs_m
  use lidort_geometry_m
  use lidort_miscsetups_m
  use lidort_la_miscsetups_m
  use lidort_lc_miscsetups_m
  use lidort_bvproblem_m
  use lidort_converge_m
  use lidort_lcs_converge_m
  use lidort_thermalsup_m
  use lidort_l_thermalsup_m
  use lidort_writemodules_m
  use lidort_l_writemodules_m
  use lidort_sfo_lcs_interface_m
  use lidort_pars_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(out) :: lidort_out_in
  type(c_ptr), intent(in) :: lidort_linfixin_in
  type(c_ptr), intent(inout) :: lidort_linmodin_in
  type(c_ptr), intent(inout) :: lidort_linsup_in
  type(c_ptr), intent(out) :: lidort_linout_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  type(lidort_fixed_lininputs), pointer :: lidort_linfixin_lcl
  type(lidort_modified_lininputs), pointer :: lidort_linmodin_lcl
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl
  type(lidort_linoutputs), pointer :: lidort_linout_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)

  call lidort_lcs_master(do_debug_input_lcl, &
                         lidort_fixin_lcl, &
                         lidort_modin_lcl, &
                         lidort_sup_lcl, &
                         lidort_out_lcl, &
                         lidort_linfixin_lcl, &
                         lidort_linmodin_lcl, &
                         lidort_linsup_lcl, &
                         lidort_linout_lcl)

end subroutine lcs_masters_m_lcs_master_wrap

end module LIDORT_LCS_MASTERS_M_WRAP

module LIDORT_LPS_MASTERS_M_WRAP

use iso_c_binding
use lidort_lps_masters_m
use lidort_pars_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine lps_masters_m_lps_master_wrap (do_debug_input_in, &
                                          lidort_fixin_in, &
                                          lidort_modin_in, &
                                          lidort_sup_in, &
                                          lidort_out_in, &
                                          lidort_linfixin_in, &
                                          lidort_linmodin_in, &
                                          lidort_linsup_in, &
                                          lidort_linout_in) bind(C)
  use lidort_pars_m
  use lidort_io_defs_m
  use lidort_lin_io_defs_m
  use lidort_inputs_m
  use lidort_l_inputs_m
  use lidort_geometry_m
  use lidort_miscsetups_m
  use lidort_la_miscsetups_m
  use lidort_lp_miscsetups_m
  use lidort_bvproblem_m
  use lidort_converge_m
  use lidort_lps_converge_m
  use lidort_thermalsup_m
  use lidort_l_thermalsup_m
  use lidort_writemodules_m
  use lidort_l_writemodules_m
  use lidort_sfo_lps_interface_m
  use lidort_pars_m

  ! Arguments
  logical(c_bool), intent(in) :: do_debug_input_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(out) :: lidort_out_in
  type(c_ptr), intent(in) :: lidort_linfixin_in
  type(c_ptr), intent(inout) :: lidort_linmodin_in
  type(c_ptr), intent(inout) :: lidort_linsup_in
  type(c_ptr), intent(out) :: lidort_linout_in

  ! Local variables
  logical(kind=4) :: do_debug_input_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  type(lidort_fixed_lininputs), pointer :: lidort_linfixin_lcl
  type(lidort_modified_lininputs), pointer :: lidort_linmodin_lcl
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl
  type(lidort_linoutputs), pointer :: lidort_linout_lcl

  ! Convert input arguments
  do_debug_input_lcl = do_debug_input_in
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)

  call lidort_lps_master(do_debug_input_lcl, &
                         lidort_fixin_lcl, &
                         lidort_modin_lcl, &
                         lidort_sup_lcl, &
                         lidort_out_lcl, &
                         lidort_linfixin_lcl, &
                         lidort_linmodin_lcl, &
                         lidort_linsup_lcl, &
                         lidort_linout_lcl)

end subroutine lps_masters_m_lps_master_wrap

end module LIDORT_LPS_MASTERS_M_WRAP

module LIDORT_BRDF_SUP_ACCESSORIES_M_WRAP

use iso_c_binding
use lidort_brdf_sup_accessories_m
use lidort_pars_m

! This module was auto-generated 

implicit none

contains

subroutine brdf_sup_accessories_m_brdf_input_check_wrap (brdf_sup_in_in, &
                                                         lidort_fixin_in, &
                                                         lidort_modin_in, &
                                                         lidort_brdfcheck_status_in) bind(C)
  use lidort_pars_m
  use brdf_sup_inputs_def_m
  use lidort_inputs_def_m
  use lidort_outputs_def_m

  ! Arguments
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(in) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_brdfcheck_status_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_exception_handling), pointer :: lidort_brdfcheck_status_lcl

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_brdfcheck_status_in, lidort_brdfcheck_status_lcl)

  call lidort_brdf_input_check(brdf_sup_in_lcl, &
                               lidort_fixin_lcl, &
                               lidort_modin_lcl, &
                               lidort_brdfcheck_status_lcl)

end subroutine brdf_sup_accessories_m_brdf_input_check_wrap

subroutine brdf_sup_accessories_m_brdf_input_check_error_wrap (errorfile_in_len, &
                                                               errorfile_in, &
                                                               lidort_brdfcheck_status_in) bind(C)
  use brdf_sup_aux_m
  use lidort_outputs_def_m

  ! Arguments
  integer(c_int), intent(in) :: errorfile_in_len
  character(kind=c_char) , intent(inout) :: errorfile_in(errorfile_in_len+1)
  type(c_ptr), intent(in) :: lidort_brdfcheck_status_in

  ! Local variables
  character(kind=c_char, len=errorfile_in_len) :: errorfile_lcl
  integer :: len_idx
  type(lidort_exception_handling), pointer :: lidort_brdfcheck_status_lcl

  ! Convert input arguments
  do len_idx = 1, errorfile_in_len
    errorfile_lcl(len_idx:len_idx) = &
      errorfile_in(len_idx)
  end do
  call c_f_pointer(lidort_brdfcheck_status_in, lidort_brdfcheck_status_lcl)

  call lidort_brdf_input_check_error(errorfile_lcl, &
                                     lidort_brdfcheck_status_lcl)

end subroutine brdf_sup_accessories_m_brdf_input_check_error_wrap

subroutine brdf_sup_accessories_m_set_brdf_inputs_wrap (brdf_sup_out_in, &
                                                        lidort_fixin_in, &
                                                        lidort_modin_in, &
                                                        lidort_sup_in) bind(C)
  use brdf_sup_outputs_def_m
  use lidort_pars_m
  use lidort_io_defs_m
  use lidort_sup_inout_def_m

  ! Arguments
  type(c_ptr), intent(in) :: brdf_sup_out_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(in) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in

  ! Local variables
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl

  ! Convert input arguments
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)

  call set_lidort_brdf_inputs(brdf_sup_out_lcl, &
                              lidort_fixin_lcl, &
                              lidort_modin_lcl, &
                              lidort_sup_lcl)

end subroutine brdf_sup_accessories_m_set_brdf_inputs_wrap

end module LIDORT_BRDF_SUP_ACCESSORIES_M_WRAP
