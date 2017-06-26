
module BRDF_LINSUP_MASTERS_M_WRAP

use iso_c_binding
use brdf_linsup_masters_m
use lidort_pars

! This module was auto-generated 

implicit none

contains

subroutine brdf_linsup_masters_m_brdf_lin_inputmaster_wrap (filnam_in_len, &
                                                            filnam_in, &
                                                            brdf_sup_in_in, &
                                                            brdf_linsup_in_in, &
                                                            brdf_sup_inputstatus_in) bind(C)
  use lidort_pars
  use brdf_findpar_m
  use brdf_sup_inputs_def
  use brdf_sup_outputs_def
  use brdf_linsup_inputs_def

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

end subroutine brdf_linsup_masters_m_brdf_lin_inputmaster_wrap

subroutine brdf_linsup_masters_m_brdf_lin_mainmaster_wrap (thread_in, &
                                                           do_debug_restoration_in, &
                                                           nmoments_input_in, &
                                                           brdf_sup_in_in, &
                                                           brdf_linsup_in_in, &
                                                           brdf_sup_out_in, &
                                                           brdf_linsup_out_in) bind(C)
  use lidort_pars
  use brdf_sup_inputs_def
  use brdf_linsup_inputs_def
  use brdf_sup_outputs_def
  use brdf_linsup_outputs_def
  use brdf_sup_aux_m
  use brdf_sup_routines_m
  use brdf_linsup_routines_m

  ! Arguments
  integer(c_int), intent(in) :: thread_in
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(in) :: brdf_linsup_in_in
  type(c_ptr), intent(out) :: brdf_sup_out_in
  type(c_ptr), intent(out) :: brdf_linsup_out_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_linsup_inputs), pointer :: brdf_linsup_in_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(brdf_linsup_outputs), pointer :: brdf_linsup_out_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_linsup_in_in, brdf_linsup_in_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(brdf_linsup_out_in, brdf_linsup_out_lcl)

  call brdf_lin_mainmaster(thread_in, &
                           do_debug_restoration_lcl, &
                           nmoments_input_in, &
                           brdf_sup_in_lcl, &
                           brdf_linsup_in_lcl, &
                           brdf_sup_out_lcl, &
                           brdf_linsup_out_lcl)

end subroutine brdf_linsup_masters_m_brdf_lin_mainmaster_wrap

end module BRDF_LINSUP_MASTERS_M_WRAP

module BRDF_SUP_MASTERS_M_WRAP

use iso_c_binding
use brdf_sup_masters_m
use lidort_pars

! This module was auto-generated 

implicit none

contains

subroutine brdf_sup_masters_m_brdf_inputmaster_wrap (filnam_in_len, &
                                                     filnam_in, &
                                                     brdf_sup_in_in, &
                                                     brdf_sup_inputstatus_in) bind(C)
  use lidort_pars
  use brdf_findpar_m
  use brdf_sup_inputs_def
  use brdf_sup_outputs_def

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

subroutine brdf_sup_masters_m_brdf_mainmaster_wrap (thread_in, &
                                                    do_debug_restoration_in, &
                                                    nmoments_input_in, &
                                                    brdf_sup_in_in, &
                                                    brdf_sup_out_in) bind(C)
  use lidort_pars
  use brdf_sup_inputs_def
  use brdf_sup_outputs_def
  use brdf_sup_aux_m
  use brdf_sup_routines_m

  ! Arguments
  integer(c_int), intent(in) :: thread_in
  logical(c_bool), intent(in) :: do_debug_restoration_in
  integer(c_int), intent(in) :: nmoments_input_in
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_sup_out_in

  ! Local variables
  logical(kind=4) :: do_debug_restoration_lcl
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl

  ! Convert input arguments
  do_debug_restoration_lcl = do_debug_restoration_in
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)

  call brdf_mainmaster(thread_in, &
                       do_debug_restoration_lcl, &
                       nmoments_input_in, &
                       brdf_sup_in_lcl, &
                       brdf_sup_out_lcl)

end subroutine brdf_sup_masters_m_brdf_mainmaster_wrap

end module BRDF_SUP_MASTERS_M_WRAP

module LIDORT_LCS_MASTERS_WRAP

use iso_c_binding
use lidort_lcs_masters
use lidort_pars
use lidort_pars
use lidort_io_defs
use lidort_lin_io_defs
use lidort_aux
use lidort_inputs
use lidort_geometry
use lidort_corrections
use lidort_lc_corrections
use lidort_ls_corrections
use lidort_miscsetups
use lidort_la_miscsetups
use lidort_lc_miscsetups
use lidort_solutions
use lidort_lpc_solutions
use lidort_bvproblem
use lidort_lc_bvproblem
use lidort_intensity
use lidort_lc_wfatmos
use lidort_ls_wfsurface
use lidort_ls_wfsleave
use lidort_thermalsup
use lidort_l_thermalsup
use lidort_writemodules
use lidort_l_writemodules

! This module was auto-generated 

implicit none

contains

subroutine lcs_masters_lcs_master_wrap (thread_in, &
                                        lidort_fixin_in, &
                                        lidort_modin_in, &
                                        lidort_sup_in, &
                                        lidort_out_in, &
                                        lidort_linfixin_in, &
                                        lidort_linmodin_in, &
                                        lidort_linsup_in, &
                                        lidort_linout_in) bind(C)
  use lidort_pars
  use lidort_pars
  use lidort_io_defs
  use lidort_lin_io_defs
  use lidort_aux
  use lidort_inputs
  use lidort_geometry
  use lidort_corrections
  use lidort_lc_corrections
  use lidort_ls_corrections
  use lidort_miscsetups
  use lidort_la_miscsetups
  use lidort_lc_miscsetups
  use lidort_solutions
  use lidort_lpc_solutions
  use lidort_bvproblem
  use lidort_lc_bvproblem
  use lidort_intensity
  use lidort_lc_wfatmos
  use lidort_ls_wfsurface
  use lidort_ls_wfsleave
  use lidort_thermalsup
  use lidort_l_thermalsup
  use lidort_writemodules
  use lidort_l_writemodules

  ! Arguments
  integer(c_int), intent(in) :: thread_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(inout) :: lidort_out_in
  type(c_ptr), intent(in) :: lidort_linfixin_in
  type(c_ptr), intent(inout) :: lidort_linmodin_in
  type(c_ptr), intent(inout) :: lidort_linsup_in
  type(c_ptr), intent(inout) :: lidort_linout_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  type(lidort_fixed_lininputs), pointer :: lidort_linfixin_lcl
  type(lidort_modified_lininputs), pointer :: lidort_linmodin_lcl
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl
  type(lidort_linoutputs), pointer :: lidort_linout_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)

  call lidort_lcs_master(thread_in, &
                         lidort_fixin_lcl, &
                         lidort_modin_lcl, &
                         lidort_sup_lcl, &
                         lidort_out_lcl, &
                         lidort_linfixin_lcl, &
                         lidort_linmodin_lcl, &
                         lidort_linsup_lcl, &
                         lidort_linout_lcl)

end subroutine lcs_masters_lcs_master_wrap

end module LIDORT_LCS_MASTERS_WRAP

module LIDORT_LPS_MASTERS_WRAP

use iso_c_binding
use lidort_lps_masters
use lidort_pars
use lidort_pars
use lidort_io_defs
use lidort_lin_io_defs
use lidort_aux
use lidort_inputs
use lidort_geometry
use lidort_corrections
use lidort_lp_corrections
use lidort_ls_corrections
use lidort_miscsetups
use lidort_la_miscsetups
use lidort_lp_miscsetups
use lidort_solutions
use lidort_lpc_solutions
use lidort_bvproblem
use lidort_lp_bvproblem
use lidort_intensity
use lidort_lp_wfatmos
use lidort_ls_wfsurface
use lidort_ls_wfsleave
use lidort_thermalsup
use lidort_l_thermalsup
use lidort_writemodules
use lidort_l_writemodules

! This module was auto-generated 

implicit none

contains

subroutine lps_masters_lps_master_wrap (thread_in, &
                                        lidort_fixin_in, &
                                        lidort_modin_in, &
                                        lidort_sup_in, &
                                        lidort_out_in, &
                                        lidort_linfixin_in, &
                                        lidort_linmodin_in, &
                                        lidort_linsup_in, &
                                        lidort_linout_in) bind(C)
  use lidort_pars
  use lidort_pars
  use lidort_io_defs
  use lidort_lin_io_defs
  use lidort_aux
  use lidort_inputs
  use lidort_geometry
  use lidort_corrections
  use lidort_lp_corrections
  use lidort_ls_corrections
  use lidort_miscsetups
  use lidort_la_miscsetups
  use lidort_lp_miscsetups
  use lidort_solutions
  use lidort_lpc_solutions
  use lidort_bvproblem
  use lidort_lp_bvproblem
  use lidort_intensity
  use lidort_lp_wfatmos
  use lidort_ls_wfsurface
  use lidort_ls_wfsleave
  use lidort_thermalsup
  use lidort_l_thermalsup
  use lidort_writemodules
  use lidort_l_writemodules

  ! Arguments
  integer(c_int), intent(in) :: thread_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(inout) :: lidort_out_in
  type(c_ptr), intent(in) :: lidort_linfixin_in
  type(c_ptr), intent(inout) :: lidort_linmodin_in
  type(c_ptr), intent(inout) :: lidort_linsup_in
  type(c_ptr), intent(inout) :: lidort_linout_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  type(lidort_fixed_lininputs), pointer :: lidort_linfixin_lcl
  type(lidort_modified_lininputs), pointer :: lidort_linmodin_lcl
  type(lidort_linsup_inout), pointer :: lidort_linsup_lcl
  type(lidort_linoutputs), pointer :: lidort_linout_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)

  call lidort_lps_master(thread_in, &
                         lidort_fixin_lcl, &
                         lidort_modin_lcl, &
                         lidort_sup_lcl, &
                         lidort_out_lcl, &
                         lidort_linfixin_lcl, &
                         lidort_linmodin_lcl, &
                         lidort_linsup_lcl, &
                         lidort_linout_lcl)

end subroutine lps_masters_lps_master_wrap

end module LIDORT_LPS_MASTERS_WRAP

module LIDORT_INPUTS_WRAP

use iso_c_binding
use lidort_inputs
use lidort_pars
use lidort_aux

! This module was auto-generated 

implicit none

contains

subroutine inputs_input_master_wrap (filnam_in_len, &
                                     filnam_in, &
                                     lidort_fixin_in, &
                                     lidort_modin_in, &
                                     lidort_inputstatus_in) bind(C)
  use lidort_pars
  use lidort_inputs_def
  use lidort_outputs_def
  use lidort_aux

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

end subroutine inputs_input_master_wrap

end module LIDORT_INPUTS_WRAP

module LIDORT_MASTERS_WRAP

use iso_c_binding
use lidort_masters
use lidort_pars
use lidort_pars
use lidort_io_defs
use lidort_aux
use lidort_inputs
use lidort_geometry
use lidort_corrections
use lidort_miscsetups
use lidort_solutions
use lidort_bvproblem
use lidort_intensity
use lidort_thermalsup
use lidort_writemodules

! This module was auto-generated 

implicit none

contains

subroutine masters_master_wrap (thread_in, &
                                lidort_fixin_in, &
                                lidort_modin_in, &
                                lidort_sup_in, &
                                lidort_out_in) bind(C)
  use lidort_pars
  use lidort_pars
  use lidort_io_defs
  use lidort_aux
  use lidort_inputs
  use lidort_geometry
  use lidort_corrections
  use lidort_miscsetups
  use lidort_solutions
  use lidort_bvproblem
  use lidort_intensity
  use lidort_thermalsup
  use lidort_writemodules

  ! Arguments
  integer(c_int), intent(in) :: thread_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(inout) :: lidort_out_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)

  call lidort_master(thread_in, &
                     lidort_fixin_lcl, &
                     lidort_modin_lcl, &
                     lidort_sup_lcl, &
                     lidort_out_lcl)

end subroutine masters_master_wrap

end module LIDORT_MASTERS_WRAP

module LIDORT_SUP_ACCESSORIES_WRAP

use iso_c_binding
use lidort_sup_accessories
use lidort_pars

! This module was auto-generated 

implicit none

contains

subroutine sup_accessories_brdf_input_checker_wrap (brdf_sup_in_in, &
                                                    lidort_fixin_in, &
                                                    lidort_modin_in, &
                                                    lidort_brdfcheck_status_in) bind(C)
  use lidort_pars
  use brdf_sup_inputs_def
  use lidort_inputs_def
  use lidort_outputs_def

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

  call lidort_brdf_input_checker(brdf_sup_in_lcl, &
                                 lidort_fixin_lcl, &
                                 lidort_modin_lcl, &
                                 lidort_brdfcheck_status_lcl)

end subroutine sup_accessories_brdf_input_checker_wrap

end module LIDORT_SUP_ACCESSORIES_WRAP
