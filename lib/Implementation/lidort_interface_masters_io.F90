module BRDF_LINSUP_MASTERS_M_IO

use iso_c_binding
use lidort_interface_types_io
use lidort_pars
use brdf_findpar_m
use brdf_sup_inputs_def
use brdf_sup_outputs_def
use brdf_linsup_inputs_def
use lidort_pars
use brdf_sup_inputs_def
use brdf_linsup_inputs_def
use brdf_sup_outputs_def
use brdf_linsup_outputs_def
use brdf_sup_aux_m
use brdf_sup_routines_m
use brdf_linsup_routines_m

! This module was auto-generated 

implicit none

contains

subroutine brdf_linsup_masters_m_read_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                            brdf_linsup_in_in, &
                                            brdf_sup_inputstatus_in, &
                                            brdf_sup_out_in, &
                                            brdf_linsup_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_linsup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: brdf_sup_out_in
  type(c_ptr), intent(out) :: brdf_linsup_out_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_linsup_inputs), pointer :: brdf_linsup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(brdf_linsup_outputs), pointer :: brdf_linsup_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_linsup_in_in, brdf_linsup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(brdf_linsup_out_in, brdf_linsup_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call brdf_linsup_masters_m_read(filename_lcl, brdf_sup_in_lcl, &
                                  brdf_linsup_in_lcl, &
                                  brdf_sup_inputstatus_lcl, &
                                  brdf_sup_out_lcl, &
                                  brdf_linsup_out_lcl)

end subroutine brdf_linsup_masters_m_read_wrap

subroutine brdf_linsup_masters_m_read (filename, brdf_sup_in_in, &
                                       brdf_linsup_in_in, &
                                       brdf_sup_inputstatus_in, &
                                       brdf_sup_out_in, &
                                       brdf_linsup_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(brdf_linsup_inputs), intent(inout), pointer :: brdf_linsup_in_in
  type(brdf_input_exception_handling), intent(inout), pointer :: brdf_sup_inputstatus_in
  type(brdf_sup_outputs), intent(inout), pointer :: brdf_sup_out_in
  type(brdf_linsup_outputs), intent(inout), pointer :: brdf_linsup_out_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_read(500, brdf_sup_in_in)
  call brdf_linsup_inputs_f_read(500, brdf_linsup_in_in)
  call brdf_input_exception_handling_f_read(500, brdf_sup_inputstatus_in)
  call brdf_sup_outputs_f_read(500, brdf_sup_out_in)
  call brdf_linsup_outputs_f_read(500, brdf_linsup_out_in)
  close(500)

end subroutine brdf_linsup_masters_m_read

subroutine brdf_linsup_masters_m_write_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                             brdf_linsup_in_in, &
                                             brdf_sup_inputstatus_in, &
                                             brdf_sup_out_in, &
                                             brdf_linsup_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_linsup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: brdf_sup_out_in
  type(c_ptr), intent(out) :: brdf_linsup_out_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_linsup_inputs), pointer :: brdf_linsup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  type(brdf_linsup_outputs), pointer :: brdf_linsup_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_linsup_in_in, brdf_linsup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  call c_f_pointer(brdf_linsup_out_in, brdf_linsup_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call brdf_linsup_masters_m_write(filename_lcl, brdf_sup_in_lcl, &
                                   brdf_linsup_in_lcl, &
                                   brdf_sup_inputstatus_lcl, &
                                   brdf_sup_out_lcl, &
                                   brdf_linsup_out_lcl)

end subroutine brdf_linsup_masters_m_write_wrap

subroutine brdf_linsup_masters_m_write (filename, brdf_sup_in_in, &
                                        brdf_linsup_in_in, &
                                        brdf_sup_inputstatus_in, &
                                        brdf_sup_out_in, &
                                        brdf_linsup_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(brdf_linsup_inputs), intent(inout), pointer :: brdf_linsup_in_in
  type(brdf_input_exception_handling), intent(inout), pointer :: brdf_sup_inputstatus_in
  type(brdf_sup_outputs), intent(inout), pointer :: brdf_sup_out_in
  type(brdf_linsup_outputs), intent(inout), pointer :: brdf_linsup_out_in
  
  
  open (500, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_write(500, brdf_sup_in_in)
  call brdf_linsup_inputs_f_write(500, brdf_linsup_in_in)
  call brdf_input_exception_handling_f_write(500, brdf_sup_inputstatus_in)
  call brdf_sup_outputs_f_write(500, brdf_sup_out_in)
  call brdf_linsup_outputs_f_write(500, brdf_linsup_out_in)
  close(500)

end subroutine brdf_linsup_masters_m_write


 
end module BRDF_LINSUP_MASTERS_M_IO

module BRDF_SUP_MASTERS_M_IO

use iso_c_binding
use lidort_interface_types_io
use lidort_pars
use brdf_findpar_m
use brdf_sup_inputs_def
use brdf_sup_outputs_def
use lidort_pars
use brdf_sup_inputs_def
use brdf_sup_outputs_def
use brdf_sup_aux_m
use brdf_sup_routines_m

! This module was auto-generated 

implicit none

contains

subroutine brdf_sup_masters_m_read_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                         brdf_sup_inputstatus_in, &
                                         brdf_sup_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: brdf_sup_out_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call brdf_sup_masters_m_read(filename_lcl, brdf_sup_in_lcl, &
                               brdf_sup_inputstatus_lcl, &
                               brdf_sup_out_lcl)

end subroutine brdf_sup_masters_m_read_wrap

subroutine brdf_sup_masters_m_read (filename, brdf_sup_in_in, &
                                    brdf_sup_inputstatus_in, &
                                    brdf_sup_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(brdf_input_exception_handling), intent(inout), pointer :: brdf_sup_inputstatus_in
  type(brdf_sup_outputs), intent(inout), pointer :: brdf_sup_out_in
  
  
  open (501, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_read(501, brdf_sup_in_in)
  call brdf_input_exception_handling_f_read(501, brdf_sup_inputstatus_in)
  call brdf_sup_outputs_f_read(501, brdf_sup_out_in)
  close(501)

end subroutine brdf_sup_masters_m_read

subroutine brdf_sup_masters_m_write_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                          brdf_sup_inputstatus_in, &
                                          brdf_sup_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: brdf_sup_in_in
  type(c_ptr), intent(out) :: brdf_sup_inputstatus_in
  type(c_ptr), intent(out) :: brdf_sup_out_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(brdf_input_exception_handling), pointer :: brdf_sup_inputstatus_lcl
  type(brdf_sup_outputs), pointer :: brdf_sup_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(brdf_sup_inputstatus_in, brdf_sup_inputstatus_lcl)
  call c_f_pointer(brdf_sup_out_in, brdf_sup_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call brdf_sup_masters_m_write(filename_lcl, brdf_sup_in_lcl, &
                                brdf_sup_inputstatus_lcl, &
                                brdf_sup_out_lcl)

end subroutine brdf_sup_masters_m_write_wrap

subroutine brdf_sup_masters_m_write (filename, brdf_sup_in_in, &
                                     brdf_sup_inputstatus_in, &
                                     brdf_sup_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(brdf_input_exception_handling), intent(inout), pointer :: brdf_sup_inputstatus_in
  type(brdf_sup_outputs), intent(inout), pointer :: brdf_sup_out_in
  
  
  open (501, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_write(501, brdf_sup_in_in)
  call brdf_input_exception_handling_f_write(501, brdf_sup_inputstatus_in)
  call brdf_sup_outputs_f_write(501, brdf_sup_out_in)
  close(501)

end subroutine brdf_sup_masters_m_write


 
end module BRDF_SUP_MASTERS_M_IO

module LIDORT_LCS_MASTERS_IO

use iso_c_binding
use lidort_interface_types_io
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

subroutine lcs_masters_read_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                                  lidort_modin_in, &
                                  lidort_sup_in, &
                                  lidort_out_in, &
                                  lidort_linfixin_in, &
                                  lidort_linmodin_in, &
                                  lidort_linsup_in, &
                                  lidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
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
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_lcs_masters_read(filename_lcl, lidort_fixin_lcl, &
                               lidort_modin_lcl, &
                               lidort_sup_lcl, &
                               lidort_out_lcl, &
                               lidort_linfixin_lcl, &
                               lidort_linmodin_lcl, &
                               lidort_linsup_lcl, &
                               lidort_linout_lcl)

end subroutine lcs_masters_read_wrap

subroutine lidort_lcs_masters_read (filename, lidort_fixin_in, &
                                    lidort_modin_in, &
                                    lidort_sup_in, &
                                    lidort_out_in, &
                                    lidort_linfixin_in, &
                                    lidort_linmodin_in, &
                                    lidort_linsup_in, &
                                    lidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  type(lidort_fixed_lininputs), intent(inout), pointer :: lidort_linfixin_in
  type(lidort_modified_lininputs), intent(inout), pointer :: lidort_linmodin_in
  type(lidort_linsup_inout), intent(inout), pointer :: lidort_linsup_in
  type(lidort_linoutputs), intent(inout), pointer :: lidort_linout_in
  
  
  open (502, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_read(502, lidort_fixin_in)
  call lidort_modified_inputs_f_read(502, lidort_modin_in)
  call lidort_sup_inout_f_read(502, lidort_sup_in)
  call lidort_outputs_f_read(502, lidort_out_in)
  call lidort_fixed_lininputs_f_read(502, lidort_linfixin_in)
  call lidort_modified_lininputs_f_read(502, lidort_linmodin_in)
  call lidort_linsup_inout_f_read(502, lidort_linsup_in)
  call lidort_linoutputs_f_read(502, lidort_linout_in)
  close(502)

end subroutine lidort_lcs_masters_read

subroutine lcs_masters_write_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                                   lidort_modin_in, &
                                   lidort_sup_in, &
                                   lidort_out_in, &
                                   lidort_linfixin_in, &
                                   lidort_linmodin_in, &
                                   lidort_linsup_in, &
                                   lidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
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
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_lcs_masters_write(filename_lcl, lidort_fixin_lcl, &
                                lidort_modin_lcl, &
                                lidort_sup_lcl, &
                                lidort_out_lcl, &
                                lidort_linfixin_lcl, &
                                lidort_linmodin_lcl, &
                                lidort_linsup_lcl, &
                                lidort_linout_lcl)

end subroutine lcs_masters_write_wrap

subroutine lidort_lcs_masters_write (filename, lidort_fixin_in, &
                                     lidort_modin_in, &
                                     lidort_sup_in, &
                                     lidort_out_in, &
                                     lidort_linfixin_in, &
                                     lidort_linmodin_in, &
                                     lidort_linsup_in, &
                                     lidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  type(lidort_fixed_lininputs), intent(inout), pointer :: lidort_linfixin_in
  type(lidort_modified_lininputs), intent(inout), pointer :: lidort_linmodin_in
  type(lidort_linsup_inout), intent(inout), pointer :: lidort_linsup_in
  type(lidort_linoutputs), intent(inout), pointer :: lidort_linout_in
  
  
  open (502, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_write(502, lidort_fixin_in)
  call lidort_modified_inputs_f_write(502, lidort_modin_in)
  call lidort_sup_inout_f_write(502, lidort_sup_in)
  call lidort_outputs_f_write(502, lidort_out_in)
  call lidort_fixed_lininputs_f_write(502, lidort_linfixin_in)
  call lidort_modified_lininputs_f_write(502, lidort_linmodin_in)
  call lidort_linsup_inout_f_write(502, lidort_linsup_in)
  call lidort_linoutputs_f_write(502, lidort_linout_in)
  close(502)

end subroutine lidort_lcs_masters_write


 
end module LIDORT_LCS_MASTERS_IO

module LIDORT_LPS_MASTERS_IO

use iso_c_binding
use lidort_interface_types_io
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

subroutine lps_masters_read_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                                  lidort_modin_in, &
                                  lidort_sup_in, &
                                  lidort_out_in, &
                                  lidort_linfixin_in, &
                                  lidort_linmodin_in, &
                                  lidort_linsup_in, &
                                  lidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
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
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_lps_masters_read(filename_lcl, lidort_fixin_lcl, &
                               lidort_modin_lcl, &
                               lidort_sup_lcl, &
                               lidort_out_lcl, &
                               lidort_linfixin_lcl, &
                               lidort_linmodin_lcl, &
                               lidort_linsup_lcl, &
                               lidort_linout_lcl)

end subroutine lps_masters_read_wrap

subroutine lidort_lps_masters_read (filename, lidort_fixin_in, &
                                    lidort_modin_in, &
                                    lidort_sup_in, &
                                    lidort_out_in, &
                                    lidort_linfixin_in, &
                                    lidort_linmodin_in, &
                                    lidort_linsup_in, &
                                    lidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  type(lidort_fixed_lininputs), intent(inout), pointer :: lidort_linfixin_in
  type(lidort_modified_lininputs), intent(inout), pointer :: lidort_linmodin_in
  type(lidort_linsup_inout), intent(inout), pointer :: lidort_linsup_in
  type(lidort_linoutputs), intent(inout), pointer :: lidort_linout_in
  
  
  open (503, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_read(503, lidort_fixin_in)
  call lidort_modified_inputs_f_read(503, lidort_modin_in)
  call lidort_sup_inout_f_read(503, lidort_sup_in)
  call lidort_outputs_f_read(503, lidort_out_in)
  call lidort_fixed_lininputs_f_read(503, lidort_linfixin_in)
  call lidort_modified_lininputs_f_read(503, lidort_linmodin_in)
  call lidort_linsup_inout_f_read(503, lidort_linsup_in)
  call lidort_linoutputs_f_read(503, lidort_linout_in)
  close(503)

end subroutine lidort_lps_masters_read

subroutine lps_masters_write_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                                   lidort_modin_in, &
                                   lidort_sup_in, &
                                   lidort_out_in, &
                                   lidort_linfixin_in, &
                                   lidort_linmodin_in, &
                                   lidort_linsup_in, &
                                   lidort_linout_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
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
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  call c_f_pointer(lidort_linfixin_in, lidort_linfixin_lcl)
  call c_f_pointer(lidort_linmodin_in, lidort_linmodin_lcl)
  call c_f_pointer(lidort_linsup_in, lidort_linsup_lcl)
  call c_f_pointer(lidort_linout_in, lidort_linout_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_lps_masters_write(filename_lcl, lidort_fixin_lcl, &
                                lidort_modin_lcl, &
                                lidort_sup_lcl, &
                                lidort_out_lcl, &
                                lidort_linfixin_lcl, &
                                lidort_linmodin_lcl, &
                                lidort_linsup_lcl, &
                                lidort_linout_lcl)

end subroutine lps_masters_write_wrap

subroutine lidort_lps_masters_write (filename, lidort_fixin_in, &
                                     lidort_modin_in, &
                                     lidort_sup_in, &
                                     lidort_out_in, &
                                     lidort_linfixin_in, &
                                     lidort_linmodin_in, &
                                     lidort_linsup_in, &
                                     lidort_linout_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  type(lidort_fixed_lininputs), intent(inout), pointer :: lidort_linfixin_in
  type(lidort_modified_lininputs), intent(inout), pointer :: lidort_linmodin_in
  type(lidort_linsup_inout), intent(inout), pointer :: lidort_linsup_in
  type(lidort_linoutputs), intent(inout), pointer :: lidort_linout_in
  
  
  open (503, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_write(503, lidort_fixin_in)
  call lidort_modified_inputs_f_write(503, lidort_modin_in)
  call lidort_sup_inout_f_write(503, lidort_sup_in)
  call lidort_outputs_f_write(503, lidort_out_in)
  call lidort_fixed_lininputs_f_write(503, lidort_linfixin_in)
  call lidort_modified_lininputs_f_write(503, lidort_linmodin_in)
  call lidort_linsup_inout_f_write(503, lidort_linsup_in)
  call lidort_linoutputs_f_write(503, lidort_linout_in)
  close(503)

end subroutine lidort_lps_masters_write


 
end module LIDORT_LPS_MASTERS_IO

module LIDORT_INPUTS_IO

use iso_c_binding
use lidort_interface_types_io
use lidort_pars
use lidort_inputs_def
use lidort_outputs_def
use lidort_aux

! This module was auto-generated 

implicit none

contains

subroutine inputs_read_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                             lidort_modin_in, &
                             lidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: lidort_fixin_in
  type(c_ptr), intent(out) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_inputstatus_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_input_exception_handling), pointer :: lidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_inputstatus_in, lidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_inputs_read(filename_lcl, lidort_fixin_lcl, &
                          lidort_modin_lcl, &
                          lidort_inputstatus_lcl)

end subroutine inputs_read_wrap

subroutine lidort_inputs_read (filename, lidort_fixin_in, &
                               lidort_modin_in, &
                               lidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_input_exception_handling), intent(inout), pointer :: lidort_inputstatus_in
  
  
  open (504, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_read(504, lidort_fixin_in)
  call lidort_modified_inputs_f_read(504, lidort_modin_in)
  call lidort_input_exception_handling_f_read(504, lidort_inputstatus_in)
  close(504)

end subroutine lidort_inputs_read

subroutine inputs_write_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                              lidort_modin_in, &
                              lidort_inputstatus_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(out) :: lidort_fixin_in
  type(c_ptr), intent(out) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_inputstatus_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_input_exception_handling), pointer :: lidort_inputstatus_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_inputstatus_in, lidort_inputstatus_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_inputs_write(filename_lcl, lidort_fixin_lcl, &
                           lidort_modin_lcl, &
                           lidort_inputstatus_lcl)

end subroutine inputs_write_wrap

subroutine lidort_inputs_write (filename, lidort_fixin_in, &
                                lidort_modin_in, &
                                lidort_inputstatus_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_input_exception_handling), intent(inout), pointer :: lidort_inputstatus_in
  
  
  open (504, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_write(504, lidort_fixin_in)
  call lidort_modified_inputs_f_write(504, lidort_modin_in)
  call lidort_input_exception_handling_f_write(504, lidort_inputstatus_in)
  close(504)

end subroutine lidort_inputs_write


 
end module LIDORT_INPUTS_IO

module LIDORT_MASTERS_IO

use iso_c_binding
use lidort_interface_types_io
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

subroutine masters_read_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                              lidort_modin_in, &
                              lidort_sup_in, &
                              lidort_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(inout) :: lidort_out_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_masters_read(filename_lcl, lidort_fixin_lcl, &
                           lidort_modin_lcl, &
                           lidort_sup_lcl, &
                           lidort_out_lcl)

end subroutine masters_read_wrap

subroutine lidort_masters_read (filename, lidort_fixin_in, &
                                lidort_modin_in, &
                                lidort_sup_in, &
                                lidort_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  
  
  open (505, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_read(505, lidort_fixin_in)
  call lidort_modified_inputs_f_read(505, lidort_modin_in)
  call lidort_sup_inout_f_read(505, lidort_sup_in)
  call lidort_outputs_f_read(505, lidort_out_in)
  close(505)

end subroutine lidort_masters_read

subroutine masters_write_wrap (filename_in, filename_in_len, lidort_fixin_in, &
                               lidort_modin_in, &
                               lidort_sup_in, &
                               lidort_out_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(inout) :: lidort_modin_in
  type(c_ptr), intent(inout) :: lidort_sup_in
  type(c_ptr), intent(inout) :: lidort_out_in

  ! Local variables
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_sup_inout), pointer :: lidort_sup_lcl
  type(lidort_outputs), pointer :: lidort_out_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_sup_in, lidort_sup_lcl)
  call c_f_pointer(lidort_out_in, lidort_out_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_masters_write(filename_lcl, lidort_fixin_lcl, &
                            lidort_modin_lcl, &
                            lidort_sup_lcl, &
                            lidort_out_lcl)

end subroutine masters_write_wrap

subroutine lidort_masters_write (filename, lidort_fixin_in, &
                                 lidort_modin_in, &
                                 lidort_sup_in, &
                                 lidort_out_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_sup_inout), intent(inout), pointer :: lidort_sup_in
  type(lidort_outputs), intent(inout), pointer :: lidort_out_in
  
  
  open (505, file=filename, form="unformatted", access="sequential")
  call lidort_fixed_inputs_f_write(505, lidort_fixin_in)
  call lidort_modified_inputs_f_write(505, lidort_modin_in)
  call lidort_sup_inout_f_write(505, lidort_sup_in)
  call lidort_outputs_f_write(505, lidort_out_in)
  close(505)

end subroutine lidort_masters_write


 
end module LIDORT_MASTERS_IO

module LIDORT_SUP_ACCESSORIES_IO

use iso_c_binding
use lidort_interface_types_io
use lidort_pars
use brdf_sup_inputs_def
use lidort_inputs_def
use lidort_outputs_def

! This module was auto-generated 

implicit none

contains

subroutine sup_accessories_read_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                      lidort_fixin_in, &
                                      lidort_modin_in, &
                                      lidort_brdfcheck_status_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(in) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_brdfcheck_status_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_exception_handling), pointer :: lidort_brdfcheck_status_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_brdfcheck_status_in, lidort_brdfcheck_status_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_sup_accessories_read(filename_lcl, brdf_sup_in_lcl, &
                                   lidort_fixin_lcl, &
                                   lidort_modin_lcl, &
                                   lidort_brdfcheck_status_lcl)

end subroutine sup_accessories_read_wrap

subroutine lidort_sup_accessories_read (filename, brdf_sup_in_in, &
                                        lidort_fixin_in, &
                                        lidort_modin_in, &
                                        lidort_brdfcheck_status_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_exception_handling), intent(inout), pointer :: lidort_brdfcheck_status_in
  
  
  open (506, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_read(506, brdf_sup_in_in)
  call lidort_fixed_inputs_f_read(506, lidort_fixin_in)
  call lidort_modified_inputs_f_read(506, lidort_modin_in)
  call lidort_exception_handling_f_read(506, lidort_brdfcheck_status_in)
  close(506)

end subroutine lidort_sup_accessories_read

subroutine sup_accessories_write_wrap (filename_in, filename_in_len, brdf_sup_in_in, &
                                       lidort_fixin_in, &
                                       lidort_modin_in, &
                                       lidort_brdfcheck_status_in) bind(c)

  ! Arguments
  integer(c_int), intent(in) :: filename_in_len
  character(kind=c_char) , intent(inout) :: filename_in(filename_in_len+1)
  type(c_ptr), intent(in) :: brdf_sup_in_in
  type(c_ptr), intent(in) :: lidort_fixin_in
  type(c_ptr), intent(in) :: lidort_modin_in
  type(c_ptr), intent(out) :: lidort_brdfcheck_status_in

  ! Local variables
  type(brdf_sup_inputs), pointer :: brdf_sup_in_lcl
  type(lidort_fixed_inputs), pointer :: lidort_fixin_lcl
  type(lidort_modified_inputs), pointer :: lidort_modin_lcl
  type(lidort_exception_handling), pointer :: lidort_brdfcheck_status_lcl
  character(len=filename_in_len) :: filename_lcl
  integer :: fn_idx

  ! Convert input arguments
  call c_f_pointer(brdf_sup_in_in, brdf_sup_in_lcl)
  call c_f_pointer(lidort_fixin_in, lidort_fixin_lcl)
  call c_f_pointer(lidort_modin_in, lidort_modin_lcl)
  call c_f_pointer(lidort_brdfcheck_status_in, lidort_brdfcheck_status_lcl)
  do fn_idx = 1, filename_in_len
    filename_lcl(fn_idx:fn_idx) = filename_in(fn_idx)
  end do

  call lidort_sup_accessories_write(filename_lcl, brdf_sup_in_lcl, &
                                    lidort_fixin_lcl, &
                                    lidort_modin_lcl, &
                                    lidort_brdfcheck_status_lcl)

end subroutine sup_accessories_write_wrap

subroutine lidort_sup_accessories_write (filename, brdf_sup_in_in, &
                                         lidort_fixin_in, &
                                         lidort_modin_in, &
                                         lidort_brdfcheck_status_in) 
  ! Arguments
  character (len=*), intent(in) :: filename
  type(brdf_sup_inputs), intent(inout), pointer :: brdf_sup_in_in
  type(lidort_fixed_inputs), intent(inout), pointer :: lidort_fixin_in
  type(lidort_modified_inputs), intent(inout), pointer :: lidort_modin_in
  type(lidort_exception_handling), intent(inout), pointer :: lidort_brdfcheck_status_in
  
  
  open (506, file=filename, form="unformatted", access="sequential")
  call brdf_sup_inputs_f_write(506, brdf_sup_in_in)
  call lidort_fixed_inputs_f_write(506, lidort_fixin_in)
  call lidort_modified_inputs_f_write(506, lidort_modin_in)
  call lidort_exception_handling_f_write(506, lidort_brdfcheck_status_in)
  close(506)

end subroutine lidort_sup_accessories_write


 
end module LIDORT_SUP_ACCESSORIES_IO

