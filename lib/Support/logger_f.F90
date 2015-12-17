module logger_wrap_m
  interface
     subroutine write_to_log(level, s, s_len) bind(c, name="lg_write_log_wrap")
       use iso_c_binding
       implicit none
       integer(c_int), intent(in) :: level, s_len
       character(kind=c_char), intent(in) :: s(*)
     end subroutine write_to_log
  end interface
end module logger_wrap_m
