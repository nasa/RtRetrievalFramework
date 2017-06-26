!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING LIDORT SUBROUTINES IN THE ORDER LISTED:
!
! SECTION #1:
!
! rossthin_function 
! rossthick_function
! lisparse_function                               
! lidense_function                                
! roujean_function                                
! hapke_function                                 
! rahman_function                                 
! coxmunk_function
! rherman_function
! breon_function
! rahman_vfunction (used by rherman_function & breon_function)
!
! SECTION #2:
!
! lisparse_function_plus                               
! lidense_function_plus                           
! hapke_function_plus                          
! rahman_function_plus                                 
! coxmunk_function_plus
! rherman_function_plus
! breon_function_plus
! rahman_vfunction_plus (used by rherman_function_plus & breon_function_plus)
!
!*******************************************************************************
!*******************************************************************************

module brdf_defs

  !PROGRAMMER: ROB SPURR WITH MINOR MODS BY MATT CHRISTI AND VIJAY NATRAJ
  !DATE LAST MODIFIED: 2/2/07

  !START MODULE

  !PRIVATE DATA & PROCEDURES   
  PRIVATE
  !PUBLIC DATA & PROCEDURES 
  PUBLIC :: &
       rossthin_function,&
       rossthick_function,&
       lisparse_function,&                            
       lidense_function,&                              
       roujean_function,&                             
       hapke_function,&                             
       rahman_function,&                             
       coxmunk_function,&
       rherman_function,&
       breon_function
  PUBLIC :: &  
       lisparse_function_plus,&                               
       lidense_function_plus,&                           
       hapke_function_plus,&                         
       rahman_function_plus,&                                 
       coxmunk_function_plus,&
       rherman_function_plus,&
       breon_function_plus

  ! ###############################################################
  ! #                                                             #
  ! #                    THE LIDORT  MODEL                        #
  ! #                                                             #
  ! #      (LInearized Discrete Ordinate Radiative Transfer)      #
  ! #       --         -        -        -         -              #
  ! #                                                             #
  ! ###############################################################

  ! ###############################################################
  ! #                                                             #
  ! #  Author :   Robert. J. D. Spurr                             #
  ! #                                                             #
  ! #  Address :  Harvard-Smithsonian Center for Astrophysics     #
  ! #             60 Garden Street                                #
  ! #             Cambridge, MA 02138, USA                        #
  ! #             Tel: (617) 496 7819                             #
  ! #                                                             #
  ! #  Email :      rspurr@cfa.harvard.edu                        #
  ! #                                                             #
  ! #  Version :    2.4_f90                                       #
  ! #  Release Date   March 2001 (F77 Version)                    #
  ! #  Release Date   November 2004 (F90 Version)                 #
  ! #                                                             #
  ! ###############################################################

CONTAINS

  !******************************************************************************
  !******************************************************************************

  ! ###############################################################
  ! #                                                             #
  ! # Subroutines in this Section                                 #
  ! #                                                             #
  ! #             rossthin_function                               #
  ! #             rossthick_function                              #
  ! #             lisparse_function                               #
  ! #             lidense_function                                #
  ! #             roujean_function                                #
  ! #             hapke_function                                  #
  ! #             rahman_function                                 #
  ! #             coxmunk_function                                #
  ! #             rherman_function                                #
  ! #             breon_function                                  #
  ! #                                                             #
  ! #             rahman_vfunction                                #
  ! #               (used by rherman_function &                   #
  ! #                        breon_function)                      #
  ! #                                                             #
  ! ###############################################################

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE rossthin_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       rossthin_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: rossthin_kernel

    ! local variables

    DOUBLE PRECISION :: ds1, ds2, cksi, sksi, ksi, func, pie

    ! initialise

    rossthin_kernel = 0.0d0
    pie = 4.0d0 *DATAN(1.0d0)

    ! function

    ds1  = xi * xj
    ds2  = sxi * sxj
    cksi = ds1 + ds2 * ckphi

    IF ( cksi > 1.0d0 ) cksi = 1.0d0

    sksi = DSQRT(1.0d0-cksi*cksi)
    ksi  = DACOS(cksi)
    func = ((0.5d0*pie-ksi)*cksi + sksi)/ds1

    rossthin_kernel = func - 0.5d0*pie

    ! finish

    RETURN
  END SUBROUTINE rossthin_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE rossthick_function ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       rossthick_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: rossthick_kernel

    ! local variables

    DOUBLE PRECISION :: ds1, ds2, ds3, cksi, sksi, ksi, func, pie

    ! initialise

    rossthick_kernel = 0.0d0
    pie = 4.0d0 *DATAN(1.0d0)

    ! function

    ds1  = xi * xj
    ds2  = sxi * sxj
    ds3  = xi  + xj

    cksi = ds1 + ds2 * ckphi
    IF ( cksi > 1.0d0 ) cksi = 1.0d0

    sksi = DSQRT(1.0d0-cksi*cksi)
    ksi  = DACOS(cksi)
    func = ((pie*0.5d0-ksi)*cksi + sksi)/ds3

    rossthick_kernel = func - (pie*0.25d0)

    ! finish

    RETURN
  END SUBROUTINE rossthick_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE lisparse_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       lisparse_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: lisparse_kernel

    ! local variables

    DOUBLE PRECISION :: x_inc, x_ref, sx_inc, sx_ref, ang_p, tx
    DOUBLE PRECISION :: t_inc, t_ref, t_inc_sq, t_ref_sq
    DOUBLE PRECISION :: cksi, delta, t, cost, sint, dsq, sintcost
    DOUBLE PRECISION :: a, b, h, r, p, q, dt1, dt2, dt2sq, qr, pie

    ! initialise
    !   -- return for special case

    lisparse_kernel       = 0.0d0
    IF ( ( xi == xj ) .AND. ( ckphi == 1.0d0 ) ) RETURN
    pie = 4.0d0 *DATAN(1.0d0)

    ! function
    ! ========

    ! .. incidence

    tx       = sxj / xj
    t_inc    = pars(2) * tx
    t_inc_sq = t_inc * t_inc
    ang_p    = DATAN ( t_inc )
    x_inc    = DCOS(ang_p)
    sx_inc   = DSIN(ang_p)

    ! .. reflection

    tx       = sxi / xi
    t_ref    = pars(2) * tx
    t_ref_sq = t_ref * t_ref
    ang_p    = DATAN ( t_ref )
    x_ref    = DCOS(ang_p)
    sx_ref   = DSIN(ang_p)

    ! ksi cosine

    cksi = x_inc * x_ref + sx_inc * sx_ref * ckphi

    ! contributions p and r

    p = ( 1.0d0 + cksi ) / x_ref
    a = ( 1.0d0 / x_inc )
    b = ( 1.0d0 / x_ref )
    r = a + b

    ! evaluate cos(t)

    dt1   = t_ref_sq + t_inc_sq
    dt2   = t_inc * t_ref
    dt2sq = dt2 * dt2
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    dsq   = delta * delta
    h     = DSQRT ( dsq + skphi * skphi * dt2sq )
    cost  = pars(1) * h / r

    ! set q function

    IF ( cost > 1.0d0 ) THEN
       q = 1.0d0
    ELSE
       t        = DACOS(cost)
       sint     = DSQRT ( 1.0d0 - cost * cost )
       sintcost = sint * cost
       q = 1.0d0 -  ( ( t - sintcost ) / pie )
    END IF

    ! set the kernel
    ! --------------

    qr = q * r 
    lisparse_kernel = 0.5d0 * p - qr

    !  finish

    RETURN
  END SUBROUTINE lisparse_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE lidense_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       lidense_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: lidense_kernel

    ! local variables

    DOUBLE PRECISION :: x_inc, x_ref, sx_inc, sx_ref, ang_p, tx
    DOUBLE PRECISION :: t_inc, t_ref, t_inc_sq, t_ref_sq
    DOUBLE PRECISION :: cksi, delta, t, cost, sint, dsq, sintcost
    DOUBLE PRECISION :: a, b, h, r, p, q, dt1, dt2, dt2sq, p_qr, pie

    ! initialise
    !   -- return for special case

    lidense_kernel       = 0.0d0
    IF ( ( xi == xj ) .AND. ( ckphi == 1.0d0 ) ) RETURN
    pie = 4.0d0 *DATAN(1.0d0)

    ! function
    ! ========

    ! .. incidence

    tx       = sxj / xj
    t_inc    = pars(2) * tx
    t_inc_sq = t_inc * t_inc
    ang_p    = DATAN ( t_inc )
    x_inc    = DCOS(ang_p)
    sx_inc   = DSIN(ang_p)

    ! .. reflection

    tx       = sxi / xi
    t_ref    = pars(2) * tx
    t_ref_sq = t_ref * t_ref
    ang_p    = DATAN ( t_ref )
    x_ref    = DCOS(ang_p)
    sx_ref   = DSIN(ang_p)

    ! ksi cosine

    cksi = x_inc * x_ref + sx_inc * sx_ref * ckphi

    ! contributions p and r

    p = ( 1.0d0 + cksi ) / x_ref
    a = ( 1.0d0 / x_inc )
    b = ( 1.0d0 / x_ref )
    r = a + b

    ! evaluate cos(t)

    dt1   = t_ref_sq + t_inc_sq
    dt2   = t_inc * t_ref
    dt2sq = dt2 * dt2
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    dsq   = delta * delta
    h     = DSQRT ( dsq + skphi * skphi * dt2sq )
    cost  = pars(1) * h / r

    ! set q function

    IF ( cost > 1.0d0 ) THEN
       q = 1.0d0
    ELSE
       t        = DACOS(cost)
       sint     = DSQRT ( 1.0d0 - cost * cost )
       sintcost = sint * cost
       q = 1.0d0 -  ( ( t - sintcost ) / pie )
    END IF

    ! set the kernel
    ! --------------

    p_qr = p / q / r 
    lidense_kernel = p_qr - 2.0d0

    ! finish

    RETURN
  END SUBROUTINE lidense_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE roujean_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       roujean_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: roujean_kernel

    ! local variables

    DOUBLE PRECISION :: ds1, ds2, ds3, txj, txi, phifac, s1, s2
    DOUBLE PRECISION :: xphi_r, cxphi_r, sxphi_r, xphi_c, pie

    ! initialise

    roujean_kernel = 0.0d0
    pie = 4.0d0 *DATAN(1.0d0)

    ! function

    xphi_c = xphi
    IF ( xphi > pie )   xphi_c = 2.0d0*pie - xphi
    IF ( xphi < 0.0d0 ) xphi_c = - xphi

    IF ( sxi < 0.0d0 ) THEN
       xphi_r  = ( pie - xphi_c )
       cxphi_r = DCOS ( xphi_r )
       sxphi_r = DSIN ( xphi_r )
       txi     =  - ( sxi / xi )
    ELSE
       txi     =  ( sxi / xi )
       xphi_r  = xphi_c
       cxphi_r = DCOS ( xphi_r )
       sxphi_r = DSIN ( xphi_r )
    END IF

    txj =  ( sxj / xj )
    ds1 = 2.0d0 * txj * txi
    ds2 = txj + txi
    ds3 = txj*txj  + txi*txi
    phifac = 0.25d0 *( ( pie - xphi_r ) * cxphi_r + sxphi_r ) / pie
    s1 = phifac * ds1
    s2 = ( ds2 + DSQRT ( ds3 - ds1 * cxphi_r ) ) / pie

    roujean_kernel = s1 - s2

    !  finish

    RETURN
  END SUBROUTINE roujean_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE hapke_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       hapke_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: hapke_kernel

    ! hapke kernel function.
    !   - new version, fresh coding
    !   - old version uses disort code; for validation (not included here)

    ! input variables:

    !    xi, sxi  : cosine/sine of angle of reflection (positive)
    !    xj, sxj  : cosine/sine of angle of incidence (positive)
    !    xphi     : difference of azimuth angles of incidence and reflection
    !    pars(1)  : single scattering albedo in hapke's bdr model
    !    pars(2)  : angular width parameter of opposition effect in hapke's model
    !    pars(3)  : empirical hot spot multiplier

    ! local variables
    !    b0_empir : empirical factor to account for the finite size of
    !               particles in hapke's bdr model
    !    b_hot    : term that accounts for the opposition effect
    !               (retroreflectance, hot spot) in hapke's bdr model
    !    ctheta   : cosine of phase angle in hapke's bdr model
    !    gamma    : albedo factor in hapke's bdr model
    !    phase    : scattering phase function in hapke's bdr model
    !    theta  : phase angle (radians); the angle between incidence and
    !             reflection directions in hapke's bdr model

    ! local variables

    DOUBLE PRECISION :: ctheta, theta, phase
    DOUBLE PRECISION :: hotspot, b0_empir, help_hot, b_hot
    DOUBLE PRECISION :: ssalbedo, gamma, reflec, func
    DOUBLE PRECISION :: help_j, term_j, help_i, term_i

    ! initialise

    hapke_kernel = 0.0d0

    ! geometrical part

    ! this is the code that is in disort - not right, i think.
    !      ctheta = xi * xj + DABS(sxi) *  DABS(sxj) * ckphi

    ctheta = xi * xj + sxi * sxj * ckphi
    IF ( ctheta > 1.0d0 ) ctheta = 1.0d0
    theta  = DACOS( ctheta )
    phase  = 1.0d0 + 0.5d0 * ctheta

    ! hot spot parameterization

    hotspot  = pars(2)
    b0_empir = pars(3)
    help_hot = hotspot + DTAN ( 0.5d0 * theta )
    b_hot    = b0_empir * hotspot / help_hot

    ! albedo parameterization

    ssalbedo = pars(1)
    gamma    = DSQRT ( 1.0d0 - ssalbedo )
    help_j   = 2.0d0 * xj
    term_j   = ( 1.0d0 + help_j ) / ( 1.0d0 + help_j * gamma )
    help_i   = 2.0d0 * xi
    term_i   = ( 1.0d0 + help_i ) / ( 1.0d0 + help_i * gamma )

    ! function

    reflec       = ssalbedo * 0.25d0 / ( xi + xj )
    func         = ( 1.0d0 + b_hot ) * phase + term_j * term_i - 1.0d0
    hapke_kernel = reflec * func

    ! finish

    RETURN
  END SUBROUTINE hapke_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE rahman_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       rahman_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: rahman_kernel

    ! local variables

    DOUBLE PRECISION :: t_inc, t_ref, dt1, dt2 
    DOUBLE PRECISION :: cxi, delta, k1_sq, fact
    DOUBLE PRECISION :: geom, phase, rfac, k0, k1, k2

    ! initialise

    rahman_kernel = 0.0d0
    IF ( xi == 0.0d0 .OR. xj == 0.0d0 ) RETURN

    ! parameters

    k0 = pars(1)
    k1 = pars(2)
    k2 = pars(3)

    ! geometrical angle xi

    cxi = xi * xj + sxi * sxj * ckphi
    IF ( cxi > 1.0d0 ) cxi = 1.0d0

    ! phase function

    k1_sq = k1 * k1
    fact  = ( 1.0d0 + k1_sq + 2.0d0 * k1 * cxi ) ** 1.5d0
    phase = ( 1.0d0 - k1_sq ) / fact

    ! delta and r-factor

    t_inc = sxi / xi
    t_ref = sxj / xj
    dt1   = t_inc*t_inc + t_ref*t_ref
    dt2   = t_inc * t_ref
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    rfac = ( 1.0d0 - k0 ) / ( 1.0d0 + delta )

    ! geom factor and kernel

    geom = ( xi * xj * ( xi + xj ) ) ** ( k2 - 1.0d0)
    rahman_kernel = k0 * phase * ( 1.0d0 + rfac ) * geom

    ! finish

    RETURN
  END SUBROUTINE rahman_function

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE coxmunk_function  ( maxpars, npars, pars, &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       coxmunk_kernel )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars
    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION, INTENT(OUT) :: coxmunk_kernel

    ! critical exponent taken out

    DOUBLE PRECISION, PARAMETER :: critexp = 88.0d0

    ! local variables

    DOUBLE PRECISION :: z, z1, z2, z2_sq_m1, h1, h2, rp, rl, xmp, pie
    DOUBLE PRECISION :: a, b, ta, argument, prob, fac1, fac2, ckphi_neg
    DOUBLE PRECISION :: s1, s2, s3, xxi, xxj, t1, t2, dcot
    DOUBLE PRECISION :: shadowi, shadowr, shadow
    DOUBLE PRECISION :: derfc

    ! initialise

    !print*
    !print*,'entering coxmunk'

    coxmunk_kernel = 0.0d0
    pie = 4.0d0 * DATAN(1.0d0)

    !  Comment. 2 February 2006. Vijay Natraj
    !  We have found in comparisons with the Giss Cox-Munk code that
    !  the input COSPHI (CKPHI) here is the negative of what we actually need,
    !  so introduce local variable which takes care of this
    !  Also removed factor of PIE in the kernel denominator
    !  This makes the output exactly same as GISS model for R(1,1)

    ckphi_neg = - ckphi

    ! function

    ! ..scatter angles

    ! old   z = - xi * xj + sxi * sxj * ckphi   
    ! old   IF ( z < -1.0d0) z = -1.0d0
    ! old   z1 = DACOS(-z)
    ! old   z2 = DCOS(z1*0.5d0)

    z = xi * xj + sxi * sxj * ckphi_neg
    IF ( z > 1.0d0) z = 1.0d0
    z1 = DACOS(z)
    z2 = DCOS(z1*0.5d0)  !lambda

    ! .. fresnel coefficients

    z2_sq_m1 = z2 * z2 - 1.0d0         !lambda**2 - 1
    h1 = pars(2) * z2                  !m**2 * lambda
    h2 = DSQRT ( pars(2) + z2_sq_m1 )  !c
    rp = ( h1 - h2 ) / ( h1 + h2 ) 
    rl = ( z2 - h2 ) / ( z2 + h2 )
    xmp = 0.5d0 * ( rp*rp + rl*rl )

    ! coxmunk function

    a = 2.0d0 * z2                 
    b = ( xi + xj ) / a                  
    IF ( b > 1.0d0 ) b = 1.0d0            
    a  = 0.5d0 * pie - DASIN(b)
    ta = DTAN(a)             
    argument = ta * ta  / pars(1)
    IF ( argument < critexp ) THEN
       prob = DEXP ( - argument )
       !     fac1 = prob / pie  / pars(1)
       fac1 = prob / pars(1)     
       fac2 = 0.25d0 / xi / ( b ** 4.0d0 )
       coxmunk_kernel = xmp * fac1 * fac2 / xj
    END IF

    !  No Shadow code if not flagged
    !  Shadow parameter changed to PARS(4), 8/16/2010, V. Natraj

    IF ( DABS(pars(4)) .ge. 1.d-8 ) then

    !  Shadow code

      s1 = DSQRT(pars(1)/pie)
      s3 = 1.d0/(DSQRT(pars(1)))
      s2 = s3*s3

      xxi  = xi*xi
      dcot = xi/DSQRT(1.d0-xxi)
      t1   = DEXP(-dcot*dcot*s2)
      t2   = derfc(dcot*s3)
      shadowi = 0.5d0*(s1*t1/dcot-t2)

      xxj  = xj*xj
      dcot = xj/DSQRT(1.d0-xxj)
      t1   = DEXP(-dcot*dcot*s2)
      t2   = derfc(dcot*s3)
      shadowr = 0.5d0*(s1*t1/dcot-t2)

      shadow =1.d0/(1.d0+shadowi+shadowr)

      coxmunk_kernel = coxmunk_kernel * shadow

    ENDIF

!  8/16/2010 Lambertian component added, V. Natraj
!  PARS(3) = ALBEDO

    coxmunk_kernel = coxmunk_kernel + pars(3)

    !print*
    !print*,'exiting coxmunk'

    ! finish

    RETURN
  END SUBROUTINE coxmunk_function

  !******************************************************************************
  !******************************************************************************
  subroutine rherman_function ( MAXPARS, NPARS, PARS, &
       XJ_IN, SXJ, XI_IN, SXI, &
       XPHI_REF, CKPHI_REF, SKPHI_REF, &
       RHERMAN_KERNEL )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER, INTENT(IN)           :: MAXPARS, NPARS
    DOUBLE PRECISION, INTENT(IN)  :: PARS ( MAXPARS )
    DOUBLE PRECISION, INTENT(IN)  :: XI_IN, SXI, XJ_IN, SXJ
    DOUBLE PRECISION, INTENT(IN)  :: XPHI_REF, CKPHI_REF, SKPHI_REF

    !  Subroutine output arguments 

    DOUBLE PRECISION, INTENT(OUT) :: RHERMAN_KERNEL

    !  Local variables      

    INTEGER, PARAMETER :: N_BRDF_STOKESSQ=16
    DOUBLE PRECISION   :: RHERMAN_VKERNEL(N_BRDF_STOKESSQ)

    INTEGER          :: O1

    DOUBLE PRECISION :: RAHMAN_VKERNEL(N_BRDF_STOKESSQ)
    DOUBLE PRECISION :: HFUNCTION

    DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
    DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
    DOUBLE PRECISION :: UNIT1, UNIT2, UNIT3, FACT1, FACTOR
    DOUBLE PRECISION :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
    DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
    DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
    DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
    DOUBLE PRECISION :: E1, E2, E3, E4
    DOUBLE PRECISION :: CF11, CF12, CF21, CF22
    DOUBLE PRECISION :: AF11, AF12, AF21, AF22
    DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
    DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

    DOUBLE PRECISION :: XI, XJ

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ 
       RHERMAN_VKERNEL(O1) = 0.0D0
    ENDDO

    !#####################################################################
    !#####################################################################
    !  COXMUNK scalar stuff...................
    !  Also removed factor of PIE in the kernel denominator
    !   This makes the output exactly same as GISS model for R(1,1)
    !      CKPHI_NEG = - CKPHI
    !      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
    !      IF ( Z .GT. ONE) Z = ONE
    !      Z1 = DACOS(Z)
    !      Z2 = DCOS(Z1*HALF)
    !  .. Fresnel coefficients
    !      REFIDX = 1.5d0
    !      REFIDX_SQ = REFIDX * REFIDX
    !      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
    !      H1 = REFIDX_SQ * Z2
    !      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
    !      RP = ( H1 - H2 ) / ( H1 + H2 )
    !      RL = ( Z2 - H2 ) / ( Z2 + H2 )
    !  XMP = r(i) = (1,1) reflection coefficient for natural lights
    !      XMP = HALF * ( RP*RP + RL*RL )
    !  Comment. This is the Lenoble SOS code
    !      xind=1.50
    !      xx=cos(pi*thd/360.)
    !      yy=sqrt(xind*xind+xx*xx-1.0)
    !      zz=xx*xind*xind
    !      rl=(zz-yy)/(zz+yy)
    !      rr=(xx-yy)/(xx+yy)
    !  where
    !     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
    !     xx   = Z2
    !     yy   = H2
    !     zz   = H1
    !     xind = REFIDX
    !     rl   = RP
    !     rr   = RL
    !  Coxmunk Function
    !      A = TWO * Z2
    !      B = ( XI + XJ ) / A
    !      IF ( B .GT. ONE ) B = ONE
    !      A = PIO2 - DASIN(B)
    !      TA = DTAN(A)
    !      ARGUMENT = TA * TA  / PARS(1)
    !      IF ( ARGUMENT .LT. CRITEXP ) THEN
    !        PROB = DEXP ( - ARGUMENT )
    !        FAC1 = PROB / PARS(1)
    !        FAC2 = QUARTER / XI / ( B ** FOUR )
    !        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
    !      ENDIF
    !#####################################################################
    !#####################################################################

    !  Transcription of the RMATR subroutine from Mishchenko/Travis code.

    !   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
    !   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
    !   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
    !   CN1 AND CN2, RESPECTIVELY.

    !   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
    !   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
    !   SXI and SXJ are the respective SINES (input)
    !   XPHI_REF = REFLECTION AZIMUTH ANGLE
    !   RHERMAN_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

    !  For real case, incident azimuth taken to be zero

    XPHI_INC  = 0.0D0
    CKPHI_INC = 1.0D0
    SKPHI_INC = 0.0D0

    !  Check for limiting cases

    IF (DABS(XI_IN-1.0D0) .LT. 1.0D-9) THEN
       XI = 0.999999999999d0
    ELSE
       XI = XI_IN
    END IF

    IF (DABS(XJ_IN-1.0D0) .LT. 1.0D-9) THEN
       XJ = 0.999999999999d0
    ELSE
       XJ = XJ_IN
    END IF

    !  help variables (coordinate transformations)

    VI1 = SXI * CKPHI_INC
    VI2 = SXI * SKPHI_INC
    VI3 = -XI
    VR1 = SXJ * CKPHI_REF
    VR2 = SXJ * SKPHI_REF
    VR3 = XJ

    !    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

    UNIT1  = VI1-VR1
    UNIT2  = VI2-VR2
    UNIT3  = VI3-VR3
    FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR = DSQRT(1.0D0/FACT1)

    !   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
    !   ---------------------------------------------------------

    CN1 = 1.0D0
    CN2 = 1.5d0

    !  this is the original code, but now C-variables are real

    XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
    CXI2 = 1.0D0 - (1.0D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
    CXI2 = DSQRT(CXI2)
    C1 = CN1*XI1
    C2 = CN2*CXI2
    CRPER = (C1-C2)/(C1+C2)
    C1 = CN2*XI1
    C2 = CN1*CXI2
    CRPAR = (C1-C2)/(C1+C2)

    !  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
    !  ----------------------------------------------

    TI1 = - XI * CKPHI_INC
    TI2 = - XI * SKPHI_INC
    TI3 = - SXI

    TR1 = + XJ * CKPHI_REF
    TR2 = + XJ * SKPHI_REF
    TR3 = -SXJ

    PI1  = - SKPHI_INC
    PII2 = + CKPHI_INC
    PI3  = 0.0D0

    PR1 = - SKPHI_REF
    PR2 = + CKPHI_REF
    PR3 = 0.0D0

    PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI

    !  Set the H-function

    HFUNCTION = 0.25d0 / ( XI + XJ )

    !  Settting (1,1), (1,2), (2,1) and (2,2) components
    !  -----------------------------------------------

    CF11 =  E1*CRPER+E2*CRPAR
    CF12 = -E3*CRPER+E4*CRPAR
    CF21 = -E4*CRPER+E3*CRPAR
    CF22 =  E2*CRPER+E1*CRPAR

    AF11 = DABS(CF11)
    AF12 = DABS(CF12)
    AF21 = DABS(CF21)
    AF22 = DABS(CF22)
    AF11 = AF11*AF11
    AF12 = AF12*AF12
    AF21 = AF21*AF21
    AF22 = AF22*AF22

    FACTOR = 0.5D0 * HFUNCTION
    RHERMAN_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
    RHERMAN_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
    RHERMAN_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
    RHERMAN_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

    !  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
    !  --------------------------------------------------------------

    C21 = CF21
    C22 = CF22
    CTTTP=CF11*CF12
    CTTPT=CF11*C21
    CTTPP=CF11*C22
    CTPPT=CF12*C21
    CTPPP=CF12*C22
    CPTPP=CF21*C22

    FACTOR = HFUNCTION
    RHERMAN_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
    RHERMAN_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
    RHERMAN_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
    RHERMAN_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
    RHERMAN_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
    RHERMAN_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

    !  Add the Diffuse term to R(1,1)
    !  ------------------------------

    !   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
    !   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

    !  This is just the Rahman Kernel.........different name !!

    XPHI_RAH = 180.0D0 - XPHI_REF
    CALL rahman_vfunction ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS, &
         XJ, SXJ, XI, SXI, &
         XPHI_RAH, CKPHI_REF, SKPHI_REF, &
         RAHMAN_VKERNEL )

    !  Add to the specular term

    ! mick
    RHERMAN_KERNEL = RHERMAN_VKERNEL(1) + RAHMAN_VKERNEL(1)

    !  Here is the original code from France........
    !    cs: cosinus of the solar zenith angle
    !    cv: cosinus of the viewing zenith angle
    !    phi: azimuth between the solar and observation vertical planes
    !      xx=cs**(vk-1.)
    !      yy=cs**(vk-1.)
    !      zz=(cs+cv)**(1.-vk)
    !      FF1=rho*xx*yy/zz
    !      xx=sqrt(1-cs*cs)
    !      yy=sqrt(1-cv*cv)
    !      ww=cs*cv-xx*yy*cos(pi*phi/180.)
    !      aa=1+gt*gt+2*gt*ww
    !      FF2=(1-gt*gt)/(aa**1.5)
    !      vv=xx/cs
    !      ww=yy/cv
    !      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
    !      FF3=1+(1-rho)/(1+G)
    !      RBD=FF1*FF2*FF3
    ! RBD...........................IS THE REFLECTION COEFFICIENT

    !  Finish

    RETURN
  end subroutine rherman_function

  !******************************************************************************
  !******************************************************************************
  subroutine breon_function ( MAXPARS, NPARS, PARS, &
       XJ_IN, SXJ, XI_IN, SXI, &
       XPHI_REF, CKPHI_REF, SKPHI_REF, &
       BREON_KERNEL )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER, INTENT(IN)           :: MAXPARS, NPARS
    DOUBLE PRECISION, INTENT(IN)  :: PARS ( MAXPARS )
    DOUBLE PRECISION, INTENT(IN)  :: XI_IN, SXI, XJ_IN, SXJ
    DOUBLE PRECISION, INTENT(IN)  :: XPHI_REF, CKPHI_REF, SKPHI_REF

    !  Subroutine output arguments

    DOUBLE PRECISION, INTENT(OUT) :: BREON_KERNEL

    !  Local variables

    INTEGER, PARAMETER :: N_BRDF_STOKESSQ=16
    DOUBLE PRECISION   :: BREON_VKERNEL(N_BRDF_STOKESSQ)

    INTEGER          :: O1

    DOUBLE PRECISION :: RAHMAN_VKERNEL(N_BRDF_STOKESSQ)
    DOUBLE PRECISION :: HFUNCTION

    DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
    DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
    DOUBLE PRECISION :: UNIT1, UNIT2, UNIT3, FACT1, FACTOR
    DOUBLE PRECISION :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
    DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
    DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
    DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
    DOUBLE PRECISION :: E1, E2, E3, E4
    DOUBLE PRECISION :: CF11, CF12, CF21, CF22
    DOUBLE PRECISION :: AF11, AF12, AF21, AF22
    DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
    DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

    DOUBLE PRECISION :: XI, XJ

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ 
       BREON_VKERNEL(O1) = 0.0D0
    ENDDO

    !#####################################################################
    !#####################################################################
    !  COXMUNK scalar stuff...................
    !  Also removed factor of PIE in the kernel denominator
    !   This makes the output exactly same as GISS model for R(1,1)
    !      CKPHI_NEG = - CKPHI
    !      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
    !      IF ( Z .GT. ONE) Z = ONE
    !      Z1 = DACOS(Z)
    !      Z2 = DCOS(Z1*HALF)
    !  .. Fresnel coefficients
    !      REFIDX = 1.5d0
    !      REFIDX_SQ = REFIDX * REFIDX
    !      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
    !      H1 = REFIDX_SQ * Z2
    !      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
    !      RP = ( H1 - H2 ) / ( H1 + H2 )
    !      RL = ( Z2 - H2 ) / ( Z2 + H2 )
    !  XMP = r(i) = (1,1) reflection coefficient for natural lights
    !      XMP = HALF * ( RP*RP + RL*RL )
    !  Comment. This is the Lenoble SOS code
    !      xind=1.50
    !      xx=cos(pi*thd/360.)
    !      yy=sqrt(xind*xind+xx*xx-1.0)
    !      zz=xx*xind*xind
    !      rl=(zz-yy)/(zz+yy)
    !      rr=(xx-yy)/(xx+yy)
    !  where
    !     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
    !     xx   = Z2
    !     yy   = H2
    !     zz   = H1
    !     xind = REFIDX
    !     rl   = RP
    !     rr   = RL
    !  Coxmunk Function
    !      A = TWO * Z2
    !      B = ( XI + XJ ) / A
    !      IF ( B .GT. ONE ) B = ONE
    !      A = PIO2 - DASIN(B)
    !      TA = DTAN(A)
    !      ARGUMENT = TA * TA  / PARS(1)
    !      IF ( ARGUMENT .LT. CRITEXP ) THEN
    !        PROB = DEXP ( - ARGUMENT )
    !        FAC1 = PROB / PARS(1)
    !        FAC2 = QUARTER / XI / ( B ** FOUR )
    !        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
    !      ENDIF
    !#####################################################################
    !#####################################################################

    !  Transcription of the RMATR subroutine from Mishchenko/Travis code.

    !   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
    !   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
    !   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
    !   CN1 AND CN2, RESPECTIVELY. 
    !   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
    !   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
    !   SXI and SXJ are the respective SINES (input)
    !   XPHI_REF = REFLECTION AZIMUTH ANGLE
    !   BREON_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

    !  For real case, incident azimuth taken to be zero

    XPHI_INC  = 0.0D0
    CKPHI_INC = 1.0D0
    SKPHI_INC = 0.0D0

    !  Check for limiting cases

    IF( DABS(XI_IN-1.0D0) .LT. 1.0D-9) THEN
       XI = 0.999999999999d0
    ELSE
       XI = XI_IN
    END IF

    IF( DABS(XJ_IN-1.0D0) .LT. 1.0D-9) THEN
       XJ = 0.999999999999d0
    ELSE
       XJ = XJ_IN
    END IF

    !  help variables (coordinate transformations)

    VI1 = SXI * CKPHI_INC
    VI2 = SXI * SKPHI_INC
    VI3 = -XI
    VR1 = SXJ * CKPHI_REF
    VR2 = SXJ * SKPHI_REF
    VR3 = XJ

    !    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

    UNIT1  = VI1-VR1
    UNIT2  = VI2-VR2
    UNIT3  = VI3-VR3
    FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR = DSQRT(1.0D0/FACT1)

    !   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
    !   ---------------------------------------------------------

    CN1 = 1.0D0
    CN2 = 1.5d0

    !  this is the original code, but now C-variables are real

    XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
    CXI2 = 1.0D0 - (1.0D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
    CXI2 = DSQRT(CXI2)
    C1 = CN1*XI1
    C2 = CN2*CXI2
    CRPER = (C1-C2)/(C1+C2)
    C1 = CN2*XI1
    C2 = CN1*CXI2
    CRPAR = (C1-C2)/(C1+C2)

    !  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
    !  ----------------------------------------------

    TI1 = - XI * CKPHI_INC
    TI2 = - XI * SKPHI_INC
    TI3 = - SXI

    TR1 = + XJ * CKPHI_REF
    TR2 = + XJ * SKPHI_REF
    TR3 = -SXJ

    PI1  = - SKPHI_INC
    PII2 = + CKPHI_INC
    PI3  = 0.0D0

    PR1 = - SKPHI_REF
    PR2 = + CKPHI_REF
    PR3 = 0.0D0

    PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI

    !  Set the H-function

    HFUNCTION = 0.25d0 / ( XI * XJ )

    !  Settting (1,1), (1,2), (2,1) and (2,2) components
    !  -----------------------------------------------

    CF11 =  E1*CRPER+E2*CRPAR
    CF12 = -E3*CRPER+E4*CRPAR
    CF21 = -E4*CRPER+E3*CRPAR
    CF22 =  E2*CRPER+E1*CRPAR

    AF11 = DABS(CF11)
    AF12 = DABS(CF12)
    AF21 = DABS(CF21)
    AF22 = DABS(CF22)
    AF11 = AF11*AF11
    AF12 = AF12*AF12
    AF21 = AF21*AF21
    AF22 = AF22*AF22

    FACTOR = 0.5D0 * HFUNCTION
    BREON_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
    BREON_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
    BREON_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
    BREON_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

    !  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
    !  --------------------------------------------------------------

    C21 = CF21
    C22 = CF22
    CTTTP=CF11*CF12
    CTTPT=CF11*C21
    CTTPP=CF11*C22
    CTPPT=CF12*C21
    CTPPP=CF12*C22
    CPTPP=CF21*C22

    FACTOR = HFUNCTION
    BREON_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
    BREON_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
    BREON_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
    BREON_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
    BREON_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
    BREON_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

    !  Add the Diffuse term to R(1,1)
    !  ------------------------------

    !   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
    !   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

    !  This is just the Rahman Kernel.........different name !!

    XPHI_RAH = 180.0d0 - XPHI_REF
    CALL rahman_vfunction ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS, &
         XJ, SXJ, XI, SXI, &
         XPHI_RAH, CKPHI_REF, SKPHI_REF, &
         RAHMAN_VKERNEL )

    !  Add to the specular term

    ! mick
    BREON_KERNEL = BREON_VKERNEL(1) + RAHMAN_VKERNEL(1)

    !  Here is the original code from France........
    !    cs: cosinus of the solar zenith angle
    !    cv: cosinus of the viewing zenith angle
    !    phi: azimuth between the solar and observation vertical planes
    !      xx=cs**(vk-1.)
    !      yy=cs**(vk-1.)
    !      zz=(cs+cv)**(1.-vk)
    !      FF1=rho*xx*yy/zz
    !      xx=sqrt(1-cs*cs)
    !      yy=sqrt(1-cv*cv)
    !      ww=cs*cv-xx*yy*cos(pi*phi/180.)
    !      aa=1+gt*gt+2*gt*ww
    !      FF2=(1-gt*gt)/(aa**1.5)
    !      vv=xx/cs
    !      ww=yy/cv
    !      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
    !      FF3=1+(1-rho)/(1+G)
    !      RBD=FF1*FF2*FF3
    ! RBD...........................IS THE REFLECTION COEFFICIENT

    !  Finish

    RETURN
  end subroutine breon_function

  !******************************************************************************
  !******************************************************************************
  subroutine rahman_vfunction ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS, &
       XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
       RAHMAN_VKERNEL )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER, INTENT(IN)           :: MAXPARS, NPARS, N_BRDF_STOKESSQ
    DOUBLE PRECISION, INTENT(IN)  :: PARS ( MAXPARS )
    DOUBLE PRECISION, INTENT(IN)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

    !  Subroutine output arguments      

    DOUBLE PRECISION, INTENT(OUT) :: RAHMAN_VKERNEL(N_BRDF_STOKESSQ)

    !  Local variables

    INTEGER          :: O1
    DOUBLE PRECISION :: T_INC, T_REF, DT1, DT2 
    DOUBLE PRECISION :: CXI, DELTA, K1_SQ, FACT
    DOUBLE PRECISION :: GEOM, PHASE, RFAC, K0, K1, K2
    DOUBLE PRECISION :: XPHI, CKPHI

    DOUBLE PRECISION, PARAMETER :: PIE=3.1415926535897932D0

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ 
       RAHMAN_VKERNEL(O1) = 0.0D0
    ENDDO
    XPHI  = PIE - PHI
    CKPHI = - CPHI

    !  (1,1) Function (scalar form)

    IF ( XI .EQ. 0.0D0 .OR. XJ .EQ. 0.0D0 ) RETURN

    !  parameters

    K0 = PARS(1)
    K1 = PARS(2)
    K2 = PARS(3)

    !  geometrical angle xi

    CXI = XI * XJ + SXI * SXJ * CKPHI
    IF ( CXI .GT. 1.0D0 ) CXI = 1.0D0

    !  Phase function

    K1_SQ = K1 * K1
    FACT  = ( 1.0D0 + K1_SQ + 2.0D0 * K1 * CXI ) ** 1.5D0
    PHASE = ( 1.0D0 - K1_SQ ) / FACT

    !  Delta and R-factor

    T_INC = SXI / XI
    T_REF = SXJ / XJ
    DT1   = T_INC*T_INC + T_REF*T_REF
    DT2   = T_INC * T_REF
    DELTA = DSQRT ( DT1 - 2.0D0 * DT2 * CKPHI )
    RFAC = ( 1.0D0 - K0 ) / ( 1.0D0 + DELTA )

    !  Geom factor and kernel

    GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - 1.0D0)
    RAHMAN_VKERNEL(1) = K0 * PHASE * ( 1.0D0 + RFAC ) * GEOM

    !  Other functions

    !      Placeholder

    !  Finish

    RETURN
  end subroutine rahman_vfunction

  !******************************************************************************
  !******************************************************************************

  ! ###############################################################
  ! #                                                             #
  ! # Subroutines in this Section                                 #
  ! #                                                             #
  ! #             coxmunk_function_plus                           #
  ! #             hapke_function_plus                             #
  ! #             lisparse_function_plus                          #
  ! #             lidense_function_plus                           #
  ! #             rahman_function_plus                            #
  ! #             rherman_function_plus                           #
  ! #             breon_function_plus                             #
  ! #                                                             #
  ! #             rahman_vfunction_plus                           #
  ! #               (used by rherman_function_plus &              #
  ! #                        breon_function_plus)                 #
  ! #                                                             #
  ! #                                                             #
  ! ###############################################################

  !*******************************************************************************
  !*******************************************************************************
  SUBROUTINE lisparse_function_plus &
       ( maxpars, npars, pars, do_deriv_pars,   &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       lisparse_kernel, lisparse_derivatives )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars

    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    LOGICAL,          INTENT(IN) :: do_deriv_pars ( maxpars )

    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: lisparse_kernel
    DOUBLE PRECISION,INTENT(OUT) :: lisparse_derivatives ( maxpars )

    ! local variables

    INTEGER :: j
    DOUBLE PRECISION :: x_inc, x_ref, sx_inc, sx_ref, ang_p, tx
    DOUBLE PRECISION :: t_inc, t_ref, t_inc_sq, t_ref_sq
    DOUBLE PRECISION :: cksi, delta, t, cost, sint, dsq, sintcost
    DOUBLE PRECISION :: a, b, h, r, p, q, dt1, dt2, dt2sq, qr, pie
    DOUBLE PRECISION :: a2, r2, dx_h, dx_r, dx_q, dx_p, dx_qr, dy_q

    ! initialise
    !   -- return for special case

    lisparse_kernel       = 0.0d0
    DO j = 1, npars
       lisparse_derivatives(j) = 0.0d0
    END DO
    IF ( ( xi == xj ) .AND. ( ckphi == 1.0d0 ) ) RETURN
    pie = 4.0d0 *DATAN(1.0d0)

    ! function
    ! ========

    ! .. incidence

    tx       = sxj / xj
    t_inc    = pars(2) * tx
    t_inc_sq = t_inc * t_inc
    ang_p    = DATAN ( t_inc )
    x_inc    = DCOS(ang_p)
    sx_inc   = DSIN(ang_p)

    ! .. reflection

    tx       = sxi / xi
    t_ref    = pars(2) * tx
    t_ref_sq = t_ref * t_ref
    ang_p    = DATAN ( t_ref )
    x_ref    = DCOS(ang_p)
    sx_ref   = DSIN(ang_p)

    ! ksi cosine

    cksi = x_inc * x_ref + sx_inc * sx_ref * ckphi

    ! contributions p and r

    p = ( 1.0d0 + cksi ) / x_ref
    a = ( 1.0d0 / x_inc )
    b = ( 1.0d0 / x_ref )
    r = a + b

    IF ( do_deriv_pars(2) ) THEN
       dx_r = r * ( 1.0d0 - ( 1.0d0 / a / b ) ) / pars(2)
       r2   = r * r
       a2   = a * a
       dx_p = ( p * ( 1.0d0 + a2 ) - ( r2 / b ) ) / pars(2) / a2
    END IF

    ! evaluate cos(t)

    dt1   = t_ref_sq + t_inc_sq
    dt2   = t_inc * t_ref
    dt2sq = dt2 * dt2
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    dsq   = delta * delta
    h     = DSQRT ( dsq + skphi * skphi * dt2sq )
    cost  = pars(1) * h / r

    ! set q function and its derivatives if flagged

    IF ( cost > 1.0d0 ) THEN
       q = 1.0d0
       IF ( do_deriv_pars(2) ) dx_q = 0.0d0
       IF ( do_deriv_pars(1) ) dy_q = 0.0d0
    ELSE
       t        = DACOS(cost)
       sint     = DSQRT ( 1.0d0 - cost * cost )
       sintcost = sint * cost
       q = 1.0d0 -  ( ( t - sintcost ) / pie )
       IF ( do_deriv_pars(2) ) THEN
          dx_h = ( 2.0d0 * h * h - dsq ) / h / pars(2)
          dx_q = 2.0d0 * sintcost * ( (dx_h/h) - (dx_r/r) ) / pie
       END IF
       IF ( do_deriv_pars(1) ) THEN
          dy_q = 2.0d0 * sintcost / pie / pars(1)
       END IF
    END IF

    ! set the kernel
    ! --------------

    qr = q * r 
    lisparse_kernel = 0.5d0 * p - qr

    ! Set derivatives
    ! ---------------

    IF ( do_deriv_pars(2) ) THEN
       dx_qr = dx_r * q + dx_q * r
       lisparse_derivatives(2) = 0.5d0 * dx_p - dx_qr
    END IF
    IF ( do_deriv_pars(1) ) THEN
       lisparse_derivatives(1) = - r * dy_q
    END IF

    !  finish

    RETURN
  END SUBROUTINE lisparse_function_plus

  !*******************************************************************************
  !*******************************************************************************
  SUBROUTINE lidense_function_plus &
       ( maxpars, npars, pars, do_deriv_pars,   &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       lidense_kernel, lidense_derivatives )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars

    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    LOGICAL,          INTENT(IN) :: do_deriv_pars ( maxpars )

    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: lidense_kernel
    DOUBLE PRECISION,INTENT(OUT) :: lidense_derivatives ( maxpars )

    ! local variables

    INTEGER          :: j
    DOUBLE PRECISION :: x_inc, x_ref, sx_inc, sx_ref, ang_p, tx
    DOUBLE PRECISION :: t_inc, t_ref, t_inc_sq, t_ref_sq
    DOUBLE PRECISION :: cksi, delta, t, cost, sint, dsq, sintcost
    DOUBLE PRECISION :: a, b, h, r, p, q, dt1, dt2, dt2sq, p_qr, pie
    DOUBLE PRECISION :: a2, r2, dx_h, dx_r, dx_q, dx_p, dx_p_qr, dy_q

    ! initialise
    !   -- return for special case

    lidense_kernel       = 0.0d0
    DO j = 1, npars
       lidense_derivatives(j) = 0.0d0
    END DO
    IF ( ( xi == xj ) .AND. ( ckphi == 1.0d0 ) ) RETURN
    pie = 4.0d0 *DATAN(1.0d0)

    ! function
    ! ========

    ! .. incidence

    tx       = sxj / xj
    t_inc    = pars(2) * tx
    t_inc_sq = t_inc * t_inc
    ang_p    = DATAN ( t_inc )
    x_inc    = DCOS(ang_p)
    sx_inc   = DSIN(ang_p)

    ! .. reflection

    tx       = sxi / xi
    t_ref    = pars(2) * tx
    t_ref_sq = t_ref * t_ref
    ang_p    = DATAN ( t_ref )
    x_ref    = DCOS(ang_p)
    sx_ref   = DSIN(ang_p)

    ! ksi cosine

    cksi = x_inc * x_ref + sx_inc * sx_ref * ckphi

    ! contributions p and r

    p = ( 1.0d0 + cksi ) / x_ref
    a = ( 1.0d0 / x_inc )
    b = ( 1.0d0 / x_ref )
    r = a + b

    IF ( do_deriv_pars(2) ) THEN
       dx_r = r * ( 1.0d0 - ( 1.0d0 / a / b ) ) / pars(2)
       r2   = r * r
       a2   = a * a
       dx_p = ( p * ( 1.0d0 + a2 ) - ( r2 / b ) ) / pars(2) / a2
    END IF

    ! evaluate cos(t)

    dt1   = t_ref_sq + t_inc_sq
    dt2   = t_inc * t_ref
    dt2sq = dt2 * dt2
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    dsq   = delta * delta
    h     = DSQRT ( dsq + skphi * skphi * dt2sq )
    cost  = pars(1) * h / r

    ! set q function and derivatives if flagged

    IF ( cost > 1.0d0 ) THEN
       q = 1.0d0
       IF ( do_deriv_pars(2) ) dx_q = 0.0d0
       IF ( do_deriv_pars(1) ) dy_q = 0.0d0
    ELSE
       t        = DACOS(cost)
       sint     = DSQRT ( 1.0d0 - cost * cost )
       sintcost = sint * cost
       q = 1.0d0 -  ( ( t - sintcost ) / pie )
       IF ( do_deriv_pars(2) ) THEN
          dx_h = ( 2.0d0 * h * h - dsq ) / h / pars(2)
          dx_q = 2.0d0 * sintcost * ( (dx_h/h) - (dx_r/r) ) / pie
       END IF
       IF ( do_deriv_pars(1) ) THEN
          dy_q = 2.0d0 * sintcost / pie / pars(1)
       END IF
    END IF

    ! set the kernel
    ! --------------

    p_qr = p / q / r 
    lidense_kernel = p_qr - 2.0d0

    ! Set derivatives
    ! ---------------

    IF ( do_deriv_pars(1) ) THEN
       lidense_derivatives(1) = - p_qr * dy_q / q 
    END IF

    IF ( do_deriv_pars(2) ) THEN
       dx_p_qr = ( dx_p / p ) - ( dx_r / r ) - ( dx_q / q )
       lidense_derivatives(2) = p_qr * dx_p_qr
    END IF

    ! finish

    RETURN
  END SUBROUTINE lidense_function_plus

  !*******************************************************************************
  !*******************************************************************************
  SUBROUTINE hapke_function_plus &
       ( maxpars, npars, pars, do_deriv_pars,   &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       hapke_kernel, hapke_derivatives )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars

    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    LOGICAL,          INTENT(IN) :: do_deriv_pars ( maxpars )

    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: hapke_kernel
    DOUBLE PRECISION,INTENT(OUT) :: hapke_derivatives ( maxpars )

    ! hapke kernel function.
    !   - new version, fresh coding
    !   - old version uses disort code; for validation (not included here)

    ! input variables:

    !    xi, sxi  : cosine/sine of angle of reflection (positive)
    !    xj, sxj  : cosine/sine of angle of incidence (positive)
    !    xphi     : difference of azimuth angles of incidence and reflection
    !    pars(1)  : single scattering albedo in hapke's bdr model
    !    pars(2)  : angular width parameter of opposition effect in hapke's model
    !    pars(3)  : empirical hot spot multiplier

    ! local variables
    !    b0_empir : empirical factor to account for the finite size of
    !               particles in hapke's bdr model
    !    b_hot    : term that accounts for the opposition effect
    !               (retroreflectance, hot spot) in hapke's bdr model
    !    ctheta   : cosine of phase angle in hapke's bdr model
    !    gamma    : albedo factor in hapke's bdr model
    !    phase    : scattering phase function in hapke's bdr model
    !    theta  : phase angle (radians); the angle between incidence and
    !             reflection directions in hapke's bdr model

    ! local variables

    INTEGER          :: j
    DOUBLE PRECISION :: ctheta, theta, phase
    DOUBLE PRECISION :: hotspot, b0_empir, help_hot, b_hot
    DOUBLE PRECISION :: ssalbedo, gamma, reflec, func
    DOUBLE PRECISION :: help_j, term_j, ghelp_j
    DOUBLE PRECISION :: help_i, term_i, ghelp_i, ti_tj, dt1, dt2

    ! Initialise

    hapke_kernel = 0.0d0
    DO j = 1, npars
       hapke_derivatives(j) = 0.0d0
    END DO

    ! this is the code that is in disort - not right, i think.
    !      ctheta = xi * xj + DABS(sxi) *  DABS(sxj) * ckphi

    ctheta = xi * xj + sxi * sxj * ckphi
    IF ( ctheta > 1.0d0 ) ctheta = 1.0d0
    theta  = DACOS( ctheta )
    phase  = 1.0d0 + 0.5d0 * ctheta

    ! hot spot parameterization

    hotspot  = pars(2)
    b0_empir = pars(3)
    help_hot = hotspot + DTAN ( 0.5d0 * theta )
    b_hot    = b0_empir * hotspot / help_hot

    ! albedo parameterization

    ssalbedo = pars(1)
    gamma    = DSQRT ( 1.0d0 - ssalbedo )
    help_j   = 2.0d0 * xj
    ghelp_j  = ( 1.0d0 + help_j * gamma )
    term_j   = ( 1.0d0 + help_j ) / ( 1.0d0 + help_j * gamma )
    help_i   = 2.0d0 * xi
    ghelp_i  = ( 1.0d0 + help_i * gamma )
    term_i   = ( 1.0d0 + help_i ) / ( 1.0d0 + help_i * gamma )
    ti_tj    = term_j * term_i

    ! function

    reflec       = ssalbedo * 0.25d0 / ( xi + xj )
    func         = ( 1.0d0 + b_hot ) * phase + term_j * term_i - 1.0d0
    hapke_kernel = reflec * func

    !  ssalbedo derivative

    IF ( do_deriv_pars(1) ) THEN
       dt1 = hapke_kernel / ssalbedo
       dt2 = ( help_j / ghelp_j ) + ( help_i / ghelp_i )
       dt2 = dt2 * ti_tj * 0.5d0 / gamma
       hapke_derivatives(1) = dt1 + dt2 * reflec
    END IF

    !  Hotspot  derivative

    IF ( do_deriv_pars(2) ) THEN
       dt1 = b_hot * ( b0_empir - b_hot ) / b0_empir / hotspot
       hapke_derivatives(2) = dt1 * reflec * phase
    END IF

    !  empirical factor

    IF ( do_deriv_pars(3) ) THEN
       dt1 = b_hot / b0_empir 
       hapke_derivatives(3) = dt1 * reflec * phase
    END IF


    ! finish

    RETURN
  END SUBROUTINE hapke_function_plus

  !*******************************************************************************
  !*******************************************************************************
  SUBROUTINE rahman_function_plus &
       ( maxpars, npars, pars, do_deriv_pars,   &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       rahman_kernel, rahman_derivatives )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars

    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    LOGICAL,          INTENT(IN) :: do_deriv_pars ( maxpars )

    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: rahman_kernel
    DOUBLE PRECISION,INTENT(OUT) :: rahman_derivatives ( maxpars )

    ! local variables

    INTEGER          :: j
    DOUBLE PRECISION :: t_inc, t_ref, dt1, dt2 
    DOUBLE PRECISION :: cxi, delta, k1_sq, fact
    DOUBLE PRECISION :: geom, phase, rfac, k0, k1, k2
    DOUBLE PRECISION :: rfac1, d_k0, d_k1, d_k2
    DOUBLE PRECISION :: helpm, helpr, helpg, d_helpm, d_fact

    ! initialise

    rahman_kernel = 0.0d0
    DO j = 1, npars
       rahman_derivatives(j) = 0.0d0
    END DO
    IF ( xi == 0.0d0 .OR. xj == 0.0d0 ) RETURN

    ! parameters

    k0 = pars(1)
    k1 = pars(2)
    k2 = pars(3)

    ! geometrical angle xi

    cxi = xi * xj + sxi * sxj * ckphi
    IF ( cxi > 1.0d0 ) cxi = 1.0d0

    ! Phase function

    k1_sq = k1 * k1
    helpm = 1.0d0 - k1_sq 
    fact  = 1.0d0 + k1_sq + 2.0d0 * k1 * cxi
    phase = helpm / ( fact ** 1.5d0 )

    ! Delta and R-factor

    t_inc = sxi / xi
    t_ref = sxj / xj
    dt1   = t_inc*t_inc + t_ref*t_ref
    dt2   = t_inc * t_ref
    delta = DSQRT ( dt1 - 2.0d0 * dt2 * ckphi )
    helpr = 1.0d0 / ( 1.0d0 + delta )
    rfac  = ( 1.0d0 - k0 ) * helpr
    rfac1 = 1.0d0 + rfac

    ! Geom factor and kernel

    helpg = xi * xj * ( xi + xj )
    geom  = helpg ** ( k2 - 1.0d0)
    rahman_kernel = k0 * phase * rfac1 * geom

    ! K0 derivative

    IF ( do_deriv_pars(1) ) THEN
       d_k0   = ( 1.0d0 / k0 ) - ( helpr / rfac1 )
       rahman_derivatives(1) = rahman_kernel * d_k0
    END IF

    ! Phase function derivative

    IF ( do_deriv_pars(2) ) THEN
       d_fact  =   2.0d0 * k1 + 2.0d0  * cxi
       d_helpm = - 2.0d0  * k1
       d_k1    = ( d_helpm / helpm ) - 1.5d0 * ( d_fact / fact )
       rahman_derivatives(2) = rahman_kernel * d_k1
    END IF

    ! K2 derivative

    IF ( do_deriv_pars(3) ) THEN
       d_k2 = DLOG ( helpg )
       rahman_derivatives(3) = rahman_kernel * d_k2
    END IF

    ! finish

    RETURN
  END SUBROUTINE rahman_function_plus

  !******************************************************************************
  !******************************************************************************
  SUBROUTINE coxmunk_function_plus &
       ( maxpars, npars, pars, do_deriv_pars,   &
       xj, sxj, xi, sxi, xphi, ckphi, skphi,  &
       coxmunk_kernel, coxmunk_derivatives )

    IMPLICIT NONE

    ! Subroutine input arguments

    INTEGER,          INTENT(IN) :: maxpars, npars

    DOUBLE PRECISION, INTENT(IN) :: pars ( maxpars )
    LOGICAL,          INTENT(IN) :: do_deriv_pars ( maxpars )

    DOUBLE PRECISION, INTENT(IN) :: xi, sxi, xj, sxj, xphi, ckphi, skphi

    ! Subroutine output arguments

    DOUBLE PRECISION,INTENT(OUT) :: coxmunk_kernel
    DOUBLE PRECISION,INTENT(OUT) :: coxmunk_derivatives ( maxpars )

    ! critical exponent taken out

    DOUBLE PRECISION, PARAMETER :: critexp = 88.0d0

    ! local variables

    INTEGER :: j
    DOUBLE PRECISION :: z, z1, z2, z2_sq_m1, h1, h2, rp, rl, xmp
    DOUBLE PRECISION :: a, b, ta, argument, prob, fac1, fac2, pie, ckphi_neg
    DOUBLE PRECISION :: s1, s2, s3, xxi, xxj
    DOUBLE PRECISION :: t1_i, t2_i, dcot_i
    DOUBLE PRECISION :: t1_r, t2_r, dcot_r
    DOUBLE PRECISION :: shadowi, shadowr, shadow

    DOUBLE PRECISION :: h1h2, h2z2, ta_sq, dfac2, dh1, dh2, drp, drl
    DOUBLE PRECISION :: d_s1, d_s2, d_t1, d_t2
    DOUBLE PRECISION :: d_shadowi, d_shadowr, d_shadow
    DOUBLE PRECISION :: derfc

    ! initialise

    pie = 4.0d0 *DATAN(1.0d0)
    coxmunk_kernel = 0.0d0
    DO j = 1, npars
       coxmunk_derivatives(j) = 0.0d0
    END DO

    !  Comment. 2 February 2006. Vijay Natraj
    !  We have found in comparisons with the Giss Cox-Munk code that
    !  the input COSPHI (CKPHI) here is the negative of what we actually need,
    !  so introduce local variable which takes care of this
    !  Also removed factor of PIE in the kernel denominator
    !  This makes the output exactly same as GISS model for R(1,1)

    ckphi_neg = - ckphi

    ! function

    ! ..scatter angles

    ! old   z = - xi * xj + sxi * sxj * ckphi   
    ! old   IF ( z < -1.0d0) z = -1.0d0
    ! old   z1 = DACOS(-z)
    ! old   z2 = DCOS(z1*0.5d0)

    z = xi * xj + sxi * sxj * ckphi_neg
    IF ( z > 1.0d0) z = 1.0d0
    z1 = DACOS(z)
    z2 = DCOS(z1*0.5d0)

    ! .. fresnel coefficients

    z2_sq_m1 = z2 * z2 - 1.0d0
    h1 = pars(2) * z2
    h2 = DSQRT ( pars(2) + z2_sq_m1 )
    rp = ( h1 - h2 ) / ( h1 + h2 )
    rl = ( z2 - h2 ) / ( z2 + h2 )
    xmp = 0.5d0 * ( rp*rp + rl*rl )
    h1h2 = h1 + h2
    h2z2 = z2 + h2

    ! coxmunk function

    a = 2.0d0 * z2                 
    b = ( xi + xj ) / a                  
    IF ( b > 1.0d0 ) b = 1.0d0            
    a  = 0.5d0 * pie - DASIN(b)
    ta = DTAN(a)             
    ta_sq = ta * ta            
    argument = ta_sq  / pars(1)
    IF ( argument < critexp ) THEN
       prob = DEXP ( - argument )
       !     fac1 = prob / pie  / pars(1)
       fac1 = prob / pars(1)     
       fac2 = 0.25d0 / xi / ( b ** 4.0d0 )
       coxmunk_kernel = xmp * fac1 * fac2 / xj
    END IF

    !  inverse slope-squared derivative (regular term)

    IF ( do_deriv_pars(1) ) THEN
       IF ( argument < critexp ) THEN
          dfac2 = ( pars(1) - ta_sq ) / pars(1) / pars(1)
          coxmunk_derivatives(1) = - coxmunk_kernel * dfac2
       END IF
    END IF

    !  square refractive index derivative (regular term)

    IF ( do_deriv_pars(2) ) THEN
       IF ( argument < critexp ) THEN
          dh1 = z2
          dh2 = 0.5d0 / h2
          drp = ( dh1 * ( 1.0d0 - rp ) - dh2 * ( 1.0d0 + rp ) ) / h1h2
          drl =  - dh2 * ( 1.0d0 + rl ) / h2z2
          dfac2 = ( rp*drp + rl*drl ) / xmp
          coxmunk_derivatives(2) = coxmunk_kernel * dfac2
       END IF
    END IF

    !  No Shadow code if not flagged
    !  Shadow parameter changed to PARS(4), 8/16/2010, V. Natraj

    IF ( DABS(pars(4)) .ge. 1.d-8 ) then

    !  Shadow code

      s1 = DSQRT(pars(1)/pie)
      s3 = 1.d0/(DSQRT(pars(1)))
      s2 = s3*s3

      xxi     = xi*xi
      dcot_i  = xi/DSQRT(1.d0-xxi)
      t1_i    = DEXP(-dcot_i*dcot_i*s2)
      t2_i    = derfc(dcot_i*s3)
      shadowi = 0.5d0*(s1*t1_i/dcot_i-t2_i)

      xxj     = xj*xj
      dcot_r  = xj/DSQRT(1.d0-xxj)
      t1_r    = DEXP(-dcot_r*dcot_r*s2)
      t2_r    = derfc(dcot_r*s3)
      shadowr = 0.5d0*(s1*t1_r/dcot_r-t2_r)

      shadow =1.d0/(1.d0+shadowi+shadowr)

      coxmunk_kernel = coxmunk_kernel * shadow

    !  inverse slope-squared derivative (add the shadow)

      IF ( do_deriv_pars(1) ) THEN
        d_s1 = 0.5d0 / pie / s1
        d_s2 = - s2 * s2
        d_t1 = - t1_i * dcot_i * dcot_i * d_s2
        d_t2 = s2 * s2 * dcot_i * s1 * t1_i
        d_shadowi = 0.5d0 * ( d_s1*t1_i/dcot_i + s1*d_t1/dcot_i - d_t2 )

        d_t1 = - t1_r * dcot_r * dcot_r * d_s2
        d_t2 = s2 * s2 * dcot_r * s1 * t1_r
        d_shadowr = 0.5d0 * ( d_s1*t1_r/dcot_r + s1*d_t1/dcot_r - d_t2 )
        d_shadow = - shadow * shadow * ( d_shadowi + d_shadowr )

        coxmunk_derivatives(1) = coxmunk_derivatives(1) * shadow + &
                                 coxmunk_kernel * d_shadow / shadow
      END IF

    !  square refractive index derivative (add the shadow)

      IF ( do_deriv_pars(2) ) THEN
        coxmunk_derivatives(2) = coxmunk_derivatives(2) * shadow
      END IF

    ENDIF

    !  8/16/2010 Lambertian component added, V. Natraj  
    !  PARS(3) = ALBEDO

    coxmunk_kernel = coxmunk_kernel + pars(3)

    !  8/16/2010 Lambertian albedo derivative added, V. Natraj
    !  DO_DERIV_PARS(3) = ALBEDO DERIVATIVE

    if ( do_deriv_pars(3) ) then
      coxmunk_derivatives(3) = 1.d0
    endif

    ! finish

    RETURN
  END SUBROUTINE coxmunk_function_plus

  !*******************************************************************************
  !*******************************************************************************
  subroutine rherman_function_plus &
       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
       XJ_IN, SXJ, XI_IN, SXI, &
       XPHI_REF, CKPHI_REF, SKPHI_REF, &
       RHERMAN_KERNEL, RHERMAN_DERIVATIVES )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER,          INTENT(IN) :: MAXPARS, NPARS

    DOUBLE PRECISION, INTENT(IN) :: PARS ( MAXPARS )
    LOGICAL,          INTENT(IN) :: DO_DERIV_PARS ( MAXPARS )

    DOUBLE PRECISION, INTENT(IN) :: XI_IN, SXI, XJ_IN, SXJ
    DOUBLE PRECISION, INTENT(IN) :: XPHI_REF, CKPHI_REF, SKPHI_REF

    !  Subroutine output arguments      

    DOUBLE PRECISION,INTENT(OUT) :: RHERMAN_KERNEL      
    DOUBLE PRECISION,INTENT(OUT) :: RHERMAN_DERIVATIVES ( MAXPARS )

    !  Local variables

    INTEGER, PARAMETER :: N_BRDF_STOKESSQ=16
    DOUBLE PRECISION   :: RHERMAN_VKERNEL ( N_BRDF_STOKESSQ )
    DOUBLE PRECISION   :: RHERMAN_VDERIVATIVES ( N_BRDF_STOKESSQ, MAXPARS )

    INTEGER          :: J, O1

    DOUBLE PRECISION :: RAHMAN_VKERNEL (N_BRDF_STOKESSQ)
    DOUBLE PRECISION :: HFUNCTION

    DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
    DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
    DOUBLE PRECISION :: UNIT1, UNIT2, UNIT3, FACT1, FACTOR
    DOUBLE PRECISION :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
    DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
    DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
    DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
    DOUBLE PRECISION :: E1, E2, E3, E4
    DOUBLE PRECISION :: CF11, CF12, CF21, CF22
    DOUBLE PRECISION :: AF11, AF12, AF21, AF22
    DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
    DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

    DOUBLE PRECISION :: XI, XJ

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ
       RHERMAN_VKERNEL(O1) = 0.0D0
       DO J = 1, NPARS
          RHERMAN_VDERIVATIVES(O1,J) = 0.0D0
       ENDDO
    ENDDO

    !#####################################################################
    !#####################################################################
    !  COXMUNK scalar stuff...................
    !  Also removed factor of PIE in the kernel denominator
    !   This makes the output exactly same as GISS model for R(1,1)
    !      CKPHI_NEG = - CKPHI
    !      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
    !      IF ( Z .GT. ONE) Z = ONE
    !      Z1 = DACOS(Z)
    !      Z2 = DCOS(Z1*HALF)
    !  .. Fresnel coefficients
    !      REFIDX = 1.5d0
    !      REFIDX_SQ = REFIDX * REFIDX
    !      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
    !      H1 = REFIDX_SQ * Z2
    !      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
    !      RP = ( H1 - H2 ) / ( H1 + H2 )
    !      RL = ( Z2 - H2 ) / ( Z2 + H2 )
    !  XMP = r(i) = (1,1) reflection coefficient for natural lights
    !      XMP = HALF * ( RP*RP + RL*RL )
    !  Comment. This is the Lenoble SOS code
    !      xind=1.50
    !      xx=cos(pi*thd/360.)
    !      yy=sqrt(xind*xind+xx*xx-1.0)
    !      zz=xx*xind*xind
    !      rl=(zz-yy)/(zz+yy)
    !      rr=(xx-yy)/(xx+yy)
    !  where
    !     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
    !     xx   = Z2
    !     yy   = H2
    !     zz   = H1
    !     xind = REFIDX
    !     rl   = RP
    !     rr   = RL
    !  Coxmunk Function
    !      A = TWO * Z2
    !      B = ( XI + XJ ) / A
    !      IF ( B .GT. ONE ) B = ONE
    !      A = PIO2 - DASIN(B)
    !      TA = DTAN(A)
    !      ARGUMENT = TA * TA  / PARS(1)
    !      IF ( ARGUMENT .LT. CRITEXP ) THEN
    !        PROB = DEXP ( - ARGUMENT )
    !        FAC1 = PROB / PARS(1)
    !        FAC2 = QUARTER / XI / ( B ** FOUR )
    !        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
    !      ENDIF
    !#####################################################################
    !#####################################################################

    !  Transcription of the RMATR subroutine from Mishchenko/Travis code.

    !   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
    !   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
    !   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
    !   CN1 AND CN2, RESPECTIVELY.

    !   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
    !   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
    !   SXI and SXJ are the respective SINES (input)
    !   XPHI_REF = REFLECTION AZIMUTH ANGLE
    !   RHERMAN_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

    !  For real case, incident azimuth taken to be zero

    XPHI_INC  = 0.0D0
    CKPHI_INC = 1.0D0
    SKPHI_INC = 0.0D0

    !  Check for limiting cases

    IF( DABS(XI_IN-1.0D0) .LT. 1.0d-9) THEN
       XI = 0.999999999999d0
    ELSE
       XI = XI_IN
    END IF

    IF( DABS(XJ_IN-1.0D0) .LT. 1.0d-9) THEN
       XJ = 0.999999999999d0
    ELSE
       XJ = XJ_IN
    END IF

    !  help variables (coordinate transformations)

    VI1 = SXI * CKPHI_INC
    VI2 = SXI * SKPHI_INC
    VI3 = -XI
    VR1 = SXJ * CKPHI_REF
    VR2 = SXJ * SKPHI_REF
    VR3 = XJ

    !    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

    UNIT1  = VI1-VR1
    UNIT2  = VI2-VR2
    UNIT3  = VI3-VR3
    FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR = DSQRT(1.0D0/FACT1)

    !   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
    !   ---------------------------------------------------------

    CN1 = 1.0D0
    CN2 = 1.5d0

    !  this is the original code, but now C-variables are real

    XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
    CXI2 = 1.0D0 - (1.0D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
    CXI2 = DSQRT(CXI2)
    C1 = CN1*XI1
    C2 = CN2*CXI2
    CRPER = (C1-C2)/(C1+C2)
    C1 = CN2*XI1
    C2 = CN1*CXI2
    CRPAR = (C1-C2)/(C1+C2)

    !  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
    !  ----------------------------------------------

    TI1 = - XI * CKPHI_INC
    TI2 = - XI * SKPHI_INC
    TI3 = - SXI

    TR1 = + XJ * CKPHI_REF
    TR2 = + XJ * SKPHI_REF
    TR3 = -SXJ

    PI1  = - SKPHI_INC
    PII2 = + CKPHI_INC
    PI3  = 0.0D0

    PR1 = - SKPHI_REF
    PR2 = + CKPHI_REF
    PR3 = 0.0D0

    PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI

    !  Set the H-function

    HFUNCTION = 0.25d0 / ( XI + XJ )

    !  Settting (1,1), (1,2), (2,1) and (2,2) components
    !  -----------------------------------------------

    CF11 =  E1*CRPER+E2*CRPAR
    CF12 = -E3*CRPER+E4*CRPAR
    CF21 = -E4*CRPER+E3*CRPAR
    CF22 =  E2*CRPER+E1*CRPAR

    AF11 = DABS(CF11)
    AF12 = DABS(CF12)
    AF21 = DABS(CF21)
    AF22 = DABS(CF22)
    AF11 = AF11*AF11
    AF12 = AF12*AF12
    AF21 = AF21*AF21
    AF22 = AF22*AF22

    FACTOR = 0.5D0 * HFUNCTION
    RHERMAN_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
    RHERMAN_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
    RHERMAN_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
    RHERMAN_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

    !  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
    !  --------------------------------------------------------------

    C21 = CF21
    C22 = CF22
    CTTTP=CF11*CF12
    CTTPT=CF11*C21
    CTTPP=CF11*C22
    CTPPT=CF12*C21
    CTPPP=CF12*C22
    CPTPP=CF21*C22

    FACTOR = HFUNCTION
    RHERMAN_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
    RHERMAN_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
    RHERMAN_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
    RHERMAN_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
    RHERMAN_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
    RHERMAN_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

    !  Add the Diffuse term to R(1,1)
    !  ------------------------------

    !   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
    !   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

    !  This is just the Rahman Kernel.........different name !!

    XPHI_RAH = 180.0d0 - XPHI_REF
    CALL rahman_vfunction_plus ( MAXPARS, N_BRDF_STOKESSQ, &
         NPARS, PARS, DO_DERIV_PARS, &
         XJ, SXJ, XI, SXI, &
         XPHI_RAH, CKPHI_REF, SKPHI_REF, &
         RAHMAN_VKERNEL, RHERMAN_VDERIVATIVES )

    !  Add to the specular term

    RHERMAN_KERNEL = RHERMAN_VKERNEL(1) + RAHMAN_VKERNEL(1)

    DO J = 1, NPARS
       RHERMAN_DERIVATIVES(J) = RHERMAN_VDERIVATIVES(1,J)
    ENDDO

    !  Here is the original code from France........
    !    cs: cosinus of the solar zenith angle
    !    cv: cosinus of the viewing zenith angle
    !    phi: azimuth between the solar and observation vertical planes
    !      xx=cs**(vk-1.)
    !      yy=cs**(vk-1.)
    !      zz=(cs+cv)**(1.-vk)
    !      FF1=rho*xx*yy/zz
    !      xx=sqrt(1-cs*cs)
    !      yy=sqrt(1-cv*cv)
    !      ww=cs*cv-xx*yy*cos(pi*phi/180.)
    !      aa=1+gt*gt+2*gt*ww
    !      FF2=(1-gt*gt)/(aa**1.5)
    !      vv=xx/cs
    !      ww=yy/cv
    !      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
    !      FF3=1+(1-rho)/(1+G)
    !      RBD=FF1*FF2*FF3
    ! RBD...........................IS THE REFLECTION COEFFICIENT

    !  Finish

    RETURN
  end subroutine rherman_function_plus

  !*******************************************************************************
  !*******************************************************************************
  subroutine breon_function_plus &
       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
       XJ_IN, SXJ, XI_IN, SXI, &
       XPHI_REF, CKPHI_REF, SKPHI_REF, &
       BREON_KERNEL, BREON_DERIVATIVES )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER,          INTENT(IN) :: MAXPARS, NPARS

    DOUBLE PRECISION, INTENT(IN) :: PARS ( MAXPARS )
    LOGICAL,          INTENT(IN) :: DO_DERIV_PARS ( MAXPARS )

    DOUBLE PRECISION, INTENT(IN) :: XI_IN, SXI, XJ_IN, SXJ
    DOUBLE PRECISION, INTENT(IN) :: XPHI_REF, CKPHI_REF, SKPHI_REF

    !  Subroutine input arguments      

    DOUBLE PRECISION,INTENT(OUT) :: BREON_KERNEL
    DOUBLE PRECISION,INTENT(OUT) :: BREON_DERIVATIVES ( MAXPARS )

    !  Local variables

    INTEGER, PARAMETER :: N_BRDF_STOKESSQ=16
    DOUBLE PRECISION   :: BREON_VKERNEL      ( N_BRDF_STOKESSQ )
    DOUBLE PRECISION   :: BREON_VDERIVATIVES ( N_BRDF_STOKESSQ, MAXPARS )

    INTEGER          :: O1, J

    DOUBLE PRECISION :: RAHMAN_VKERNEL(N_BRDF_STOKESSQ)
    DOUBLE PRECISION :: HFUNCTION

    DOUBLE PRECISION :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
    DOUBLE PRECISION :: VI1, VI2, VI3, VR1, VR2, VR3
    DOUBLE PRECISION :: UNIT1, UNIT2, UNIT3, FACT1, FACTOR
    DOUBLE PRECISION :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
    DOUBLE PRECISION :: TI1, TI2, TI3, TR1, TR2, TR3
    DOUBLE PRECISION :: PI1, PII2, PI3, PR1, PR2, PR3
    DOUBLE PRECISION :: PIKR, PRKI, TIKR, TRKI
    DOUBLE PRECISION :: E1, E2, E3, E4
    DOUBLE PRECISION :: CF11, CF12, CF21, CF22
    DOUBLE PRECISION :: AF11, AF12, AF21, AF22
    DOUBLE PRECISION :: C21, C22, CTTTP, CTTPT, CTTPP
    DOUBLE PRECISION :: CTPPT, CTPPP, CPTPP

    DOUBLE PRECISION :: XI, XJ

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ 
       BREON_VKERNEL(O1) = 0.0D0
       DO J = 1, NPARS
          BREON_VDERIVATIVES(O1,J) = 0.0D0
       ENDDO
    ENDDO

    !#####################################################################
    !#####################################################################
    !  COXMUNK scalar stuff...................
    !  Also removed factor of PIE in the kernel denominator
    !   This makes the output exactly same as GISS model for R(1,1)
    !      CKPHI_NEG = - CKPHI
    !      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
    !      IF ( Z .GT. ONE) Z = ONE
    !      Z1 = DACOS(Z)
    !      Z2 = DCOS(Z1*HALF)
    !  .. Fresnel coefficients
    !      REFIDX = 1.5d0
    !      REFIDX_SQ = REFIDX * REFIDX
    !      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
    !      H1 = REFIDX_SQ * Z2
    !      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
    !      RP = ( H1 - H2 ) / ( H1 + H2 )
    !      RL = ( Z2 - H2 ) / ( Z2 + H2 )
    !  XMP = r(i) = (1,1) reflection coefficient for natural lights
    !      XMP = HALF * ( RP*RP + RL*RL )
    !  Comment. This is the Lenoble SOS code
    !      xind=1.50
    !      xx=cos(pi*thd/360.)
    !      yy=sqrt(xind*xind+xx*xx-1.0)
    !      zz=xx*xind*xind
    !      rl=(zz-yy)/(zz+yy)
    !      rr=(xx-yy)/(xx+yy)
    !  where
    !     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
    !     xx   = Z2
    !     yy   = H2
    !     zz   = H1
    !     xind = REFIDX
    !     rl   = RP
    !     rr   = RL
    !  Coxmunk Function
    !      A = TWO * Z2
    !      B = ( XI + XJ ) / A
    !      IF ( B .GT. ONE ) B = ONE
    !      A = PIO2 - DASIN(B)
    !      TA = DTAN(A)
    !      ARGUMENT = TA * TA  / PARS(1)
    !      IF ( ARGUMENT .LT. CRITEXP ) THEN
    !        PROB = DEXP ( - ARGUMENT )
    !        FAC1 = PROB / PARS(1)
    !        FAC2 = QUARTER / XI / ( B ** FOUR )
    !        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
    !      ENDIF
    !#####################################################################
    !#####################################################################

    !  Transcription of the RMATR subroutine from Mishchenko/Travis code.

    !   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
    !   ILLUMINATION FROM ABOVE FOR
    !   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
    !   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
    !   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
    !   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

    !   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER) 

    !   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
    !   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
    !   SXI and SXJ are the respective SINES (input)
    !   XPHI_REF = REFLECTION AZIMUTH ANGLE
    !   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

    !  For real case, incident azimuth taken to be zero

    XPHI_INC  = 0.0D0
    CKPHI_INC = 1.0D0
    SKPHI_INC = 0.0D0

    !  Check for limiting cases

    IF (DABS(XI_IN-1.0D0) .LT. 1.0d-9) THEN
       XI = 0.999999999999d0
    ELSE
       XI = XI_IN
    END IF

    IF (DABS(XJ_IN-1.0D0) .LT. 1.0d-9) THEN
       XJ = 0.999999999999d0
    ELSE
       XJ = XJ_IN
    END IF

    !  help variables (coordinate transformations)

    VI1 = SXI * CKPHI_INC
    VI2 = SXI * SKPHI_INC
    VI3 = -XI
    VR1 = SXJ * CKPHI_REF
    VR2 = SXJ * SKPHI_REF
    VR3 = XJ

    !    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

    UNIT1  = VI1-VR1
    UNIT2  = VI2-VR2
    UNIT3  = VI3-VR3
    FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR = DSQRT(1.0D0/FACT1)

    !   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
    !   ---------------------------------------------------------

    CN1 = 1.0D0
    CN2 = 1.5d0

    !  this is the original code, but now C-variables are real

    XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
    CXI2 = 1.0D0 - (1.0D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
    CXI2 = DSQRT(CXI2)
    C1 = CN1*XI1
    C2 = CN2*CXI2
    CRPER = (C1-C2)/(C1+C2)
    C1 = CN2*XI1
    C2 = CN1*CXI2
    CRPAR = (C1-C2)/(C1+C2)

    !  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
    !  ----------------------------------------------

    TI1 = - XI * CKPHI_INC
    TI2 = - XI * SKPHI_INC
    TI3 = - SXI

    TR1 = + XJ * CKPHI_REF
    TR2 = + XJ * SKPHI_REF
    TR3 = -SXJ

    PI1  = - SKPHI_INC
    PII2 = + CKPHI_INC
    PI3  = 0.0D0

    PR1 = - SKPHI_REF
    PR2 = + CKPHI_REF
    PR3 = 0.0D0

    PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI

    !  Set the H-function

    HFUNCTION = 0.25d0 / ( XI * XJ )

    !  Settting (1,1), (1,2), (2,1) and (2,2) components
    !  -----------------------------------------------

    CF11 =  E1*CRPER+E2*CRPAR
    CF12 = -E3*CRPER+E4*CRPAR
    CF21 = -E4*CRPER+E3*CRPAR
    CF22 =  E2*CRPER+E1*CRPAR

    AF11 = DABS(CF11)
    AF12 = DABS(CF12)
    AF21 = DABS(CF21)
    AF22 = DABS(CF22)
    AF11 = AF11*AF11
    AF12 = AF12*AF12
    AF21 = AF21*AF21
    AF22 = AF22*AF22

    FACTOR = 0.5D0 * HFUNCTION
    BREON_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
    BREON_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
    BREON_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
    BREON_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

    !  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
    !  --------------------------------------------------------------

    C21 = CF21
    C22 = CF22
    CTTTP=CF11*CF12
    CTTPT=CF11*C21
    CTTPP=CF11*C22
    CTPPT=CF12*C21
    CTPPP=CF12*C22
    CPTPP=CF21*C22

    FACTOR = HFUNCTION
    BREON_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
    BREON_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
    BREON_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
    BREON_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
    BREON_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
    BREON_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

    !  Add the Diffuse term to R(1,1)
    !  ------------------------------

    !   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
    !   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

    !  This is just the Rahman Kernel.........different name !!

    XPHI_RAH = 180.0d0 - XPHI_REF
    CALL rahman_vfunction_plus ( MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS, &
         DO_DERIV_PARS, &
         XJ, SXJ, XI, SXI, &
         XPHI_RAH, CKPHI_REF, SKPHI_REF, &
         RAHMAN_VKERNEL, BREON_VDERIVATIVES )

    !  Add to the specular term

    BREON_KERNEL = BREON_VKERNEL(1) + RAHMAN_VKERNEL(1)

    DO J = 1, NPARS
       BREON_DERIVATIVES(J) = BREON_VDERIVATIVES(1,J)
    ENDDO

    !  Here is the original code from France........
    !    cs: cosinus of the solar zenith angle
    !    cv: cosinus of the viewing zenith angle
    !    phi: azimuth between the solar and observation vertical planes
    !      xx=cs**(vk-1.)
    !      yy=cs**(vk-1.)
    !      zz=(cs+cv)**(1.-vk)
    !      FF1=rho*xx*yy/zz
    !      xx=sqrt(1-cs*cs)
    !      yy=sqrt(1-cv*cv)
    !      ww=cs*cv-xx*yy*cos(pi*phi/180.)
    !      aa=1+gt*gt+2*gt*ww
    !      FF2=(1-gt*gt)/(aa**1.5)
    !      vv=xx/cs
    !      ww=yy/cv
    !      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
    !      FF3=1+(1-rho)/(1+G)
    !      RBD=FF1*FF2*FF3
    ! RBD...........................IS THE REFLECTION COEFFICIENT

    !  Finish

    RETURN
  end subroutine breon_function_plus

  !*******************************************************************************
  !*******************************************************************************
  subroutine rahman_vfunction_plus (MAXPARS, N_BRDF_STOKESSQ, NPARS, PARS,&
       DO_DERIV_PARS, &
       XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
       RAHMAN_VKERNEL, RAHMAN_VDERIVATIVES )

    IMPLICIT NONE

    !  Subroutine input arguments

    INTEGER, INTENT(IN)          :: MAXPARS, NPARS, N_BRDF_STOKESSQ
    DOUBLE PRECISION, INTENT(IN) :: PARS ( MAXPARS )
    LOGICAL, INTENT(IN)          :: DO_DERIV_PARS ( MAXPARS )
    DOUBLE PRECISION, INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

    !  Subroutine output arguments      

    DOUBLE PRECISION, INTENT(OUT) :: &
         RAHMAN_VKERNEL ( N_BRDF_STOKESSQ )
    DOUBLE PRECISION, INTENT(OUT) :: &
         RAHMAN_VDERIVATIVES ( N_BRDF_STOKESSQ, MAXPARS )

    !  Local variables

    INTEGER          :: J, O1
    DOUBLE PRECISION :: T_INC, T_REF, DT1, DT2 
    DOUBLE PRECISION :: CXI, DELTA, K1_SQ, FACT, K0, K1, K2
    DOUBLE PRECISION :: HELPM, HELPR, HELPG, D_HELPM, D_FACT
    DOUBLE PRECISION :: GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
    DOUBLE PRECISION :: XPHI, CKPHI

    DOUBLE PRECISION, PARAMETER :: PIE=3.1415926535897932D0

    !  Initialise

    DO O1 = 1, N_BRDF_STOKESSQ
       RAHMAN_VKERNEL(O1) = 0.0D0
       DO J = 1, NPARS
          RAHMAN_VDERIVATIVES(O1,J) = 0.0D0
       ENDDO
    ENDDO
    XPHI  = PIE - PHI
    CKPHI = - CPHI
    IF ( XI .EQ. 0.0D0 .OR. XJ .EQ. 0.0D0 ) RETURN

    !  parameters

    K0 = PARS(1)
    K1 = PARS(2)
    K2 = PARS(3)

    !  geometrical angle xi

    CXI = XI * XJ + SXI * SXJ * CKPHI
    IF ( CXI .GT. 1.0D0 ) CXI = 1.0D0

    !  Phase function

    K1_SQ = K1 * K1
    HELPM = 1.0D0 - K1_SQ 
    FACT  = 1.0D0 + K1_SQ + 2.0D0 * K1 * CXI
    PHASE = HELPM / ( FACT ** 1.5D0 )

    !  Delta and R-factor

    T_INC = SXI / XI
    T_REF = SXJ / XJ
    DT1   = T_INC*T_INC + T_REF*T_REF
    DT2   = T_INC * T_REF
    DELTA = DSQRT ( DT1 - 2.0D0 * DT2 * CKPHI )
    HELPR = 1.0D0 / ( 1.0D0 + DELTA )
    RFAC  = ( 1.0D0 - K0 ) * HELPR
    RFAC1 = 1.0D0 + RFAC

    !  Geom factor and kernel

    HELPG = XI * XJ * ( XI + XJ )
    GEOM  = HELPG ** ( K2 - 1.0D0)
    RAHMAN_VKERNEL(1) = K0 * PHASE * RFAC1 * GEOM

    !  Other vector functions

    !    P L A C E H O L D E R

    !  Scalar derivatives (1,1) matrix entry
    !  -------------------------------------

    !  K0 derivative

    IF ( DO_DERIV_PARS(1) ) THEN
       D_K0   = ( 1.0D0 / K0 ) - ( HELPR / RFAC1 )
       RAHMAN_VDERIVATIVES(1,1) = RAHMAN_VKERNEL(1) * D_K0
    ENDIF

    !  Phase function derivative

    IF ( DO_DERIV_PARS(2) ) THEN
       D_FACT  =   2.0D0 * K1 + 2.0D0 * CXI
       D_HELPM = - 2.0D0 * K1
       D_K1    = ( D_HELPM / HELPM ) - 1.5D0 * ( D_FACT / FACT )
       RAHMAN_VDERIVATIVES(1,2) = RAHMAN_VKERNEL(1) * D_K1
    ENDIF

    !  K2 derivative

    IF ( DO_DERIV_PARS(3) ) THEN
       D_K2 = DLOG ( HELPG )
       RAHMAN_VDERIVATIVES(1,3) = RAHMAN_VKERNEL(1) * D_K2
    ENDIF

    !  Finish

    RETURN
  end subroutine rahman_vfunction_plus

  !*******************************************************************************
  !*******************************************************************************

end module brdf_defs
