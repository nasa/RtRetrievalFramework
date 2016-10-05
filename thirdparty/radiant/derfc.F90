        function derfc(x)

        implicit none

        double precision, intent(in) :: x
        double precision :: derfc

! Returns the complementary error function erfc(x) with fractional error 
!everywhere less than 1.2 * 10^7.

        double precision :: t,z

        z=DABS(x)
        t=1.d0/(1.d0+0.5d0*z)
        derfc=t*DEXP(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
       t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+ &
       t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
        if (x.lt.0.d0) derfc=2.d0-derfc

!        write(220,*) x,z,t,derfc

        return
        end function derfc
