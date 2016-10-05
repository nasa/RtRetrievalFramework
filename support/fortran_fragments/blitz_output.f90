! Include this code where you want to use routines to write out
! Fortran arrays in row ordering with the proper syntax for
! feeding to Blitz++ using the stream operator in C++

module blitz_output

  implicit none

  interface write_blitz
     module procedure &
          write_blitz_1d, &
          write_blitz_2d, &
          write_blitz_3d, &
          write_blitz_4d, &
          write_blitz_5d, &
          write_blitz_6d
  end interface
  
contains

! 1D
subroutine write_blitz_1d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:) :: array

  integer :: idx1
  
  write(unit,'(i4)'), &
       size(array, 1)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     write(unit,1001) array(idx1)
  end do

  write(unit, '(a)') ']'

1001 format(1p10e20.11)
 
end subroutine write_blitz_1d

! 2D
subroutine write_blitz_2d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:,:) :: array

  integer :: idx1,idx2
  
  write(unit,'(i4," x ",i4)'), &
       size(array, 1), &
       size(array, 2)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     write(unit,1002)(array(idx1,idx2),&
          idx2=lbound(array,2), ubound(array,2))
  end do

  write(unit, '(a)') ']'

1002 format(1p10e20.11)
 
end subroutine write_blitz_2d

! 3D
subroutine write_blitz_3d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:,:,:) :: array

  integer :: idx1,idx2,idx3
  
  write(unit,'(i4," x ",i4," x ",i4)'), &
       size(array, 1), &
       size(array, 2), &
       size(array, 3)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     do idx2 = lbound(array,2), ubound(array,2)
        write(unit,1003)(array(idx1,idx2,idx3),&
             idx3=lbound(array,3), ubound(array,3))
     end do
  end do

  write(unit, '(a)') ']'

1003 format(1p10e20.11)
 
end subroutine write_blitz_3d

! 4D
subroutine write_blitz_4d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:,:,:,:) :: array

  integer :: idx1,idx2,idx3,idx4
  
  write(unit,'(i4," x ",i4," x ",i4," x ",i4)'), &
       size(array, 1), &
       size(array, 2), &
       size(array, 3), &
       size(array, 4)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     do idx2 = lbound(array,2), ubound(array,2)
        do idx3 = lbound(array,3), ubound(array,3)
           write(unit,1004)(array(idx1,idx2,idx3,idx4),&
                idx4=lbound(array,4), ubound(array,4))
        end do
     end do
  end do

  write(unit, '(a)') ']'

1004 format(1p10e20.11)
 
end subroutine write_blitz_4d

! 5D
subroutine write_blitz_5d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:,:,:,:,:) :: array

  integer :: idx1,idx2,idx3,idx4,idx5
  
  write(unit,'(i4," x ",i4," x ",i4," x ",i4," x ",i4)'), &
       size(array, 1), &
       size(array, 2), &
       size(array, 3), &
       size(array, 4), &
       size(array, 5)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     do idx2 = lbound(array,2), ubound(array,2)
        do idx3 = lbound(array,3), ubound(array,3)
           do idx4 = lbound(array,4), ubound(array,4)
              write(unit,1005)(array(idx1,idx2,idx3,idx4,idx5),&
                   idx5=lbound(array,5), ubound(array,5))
           end do
        end do
     end do
  end do

  write(unit, '(a)') ']'

1005 format(1p10e20.11)
 
end subroutine write_blitz_5d

! 6D
subroutine write_blitz_6d(unit, array)
  implicit none
  integer, intent(in) :: unit
  double precision, intent(in), dimension (:,:,:,:,:,:) :: array

  integer :: idx1,idx2,idx3,idx4,idx5,idx6
  
  write(unit,'(i4," x ",i4," x ",i4," x ",i4," x ",i4," x ",i4)'), &
       size(array, 1), &
       size(array, 2), &
       size(array, 3), &
       size(array, 4), &
       size(array, 5), &
       size(array, 6)
  write(unit, '(a)') '['

  do idx1 = lbound(array,1), ubound(array,1)
     do idx2 = lbound(array,2), ubound(array,2)
        do idx3 = lbound(array,3), ubound(array,3)
           do idx4 = lbound(array,4), ubound(array,4)
              do idx5 = lbound(array,5), ubound(array,5)
                 write(unit,1006)(array(idx1,idx2,idx3,idx4,idx5,idx6),&
                      idx6=lbound(array,6), ubound(array,6))
              end do
           end do
        end do
     end do
  end do

  write(unit, '(a)') ']'

1006 format(1p10e20.11)
 
end subroutine write_blitz_6d

end module

