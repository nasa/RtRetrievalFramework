
! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.6 F90                               #
! #  Release Date :   August 2012                           #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

module lidort_geometry

!  Parameter types

   USE LIDORT_PARS, only : fpk

!private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #                                                             #
! #            LIDORT_CHAPMAN   (called by LIDORT_MASTER)       #
! #              BEAM_GEOMETRY_PREPARE                          #
! #                                                             #
! #            OUTGOING_SPHERGEOM_FINE_UP                       #
! #            OUTGOING_SPHERGEOM_FINE_DN                       #
! #            MULTI_OUTGOING_ADJUSTGEOM                        #
! #            LOSONLY_OUTGOING_ADJUSTGEOM                      #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Revisions in Version 3.5                                   #
! #                                                             #
! #     07 October 2010.                                        #
! #       Original code had only one geometry routine for the   #
! #       Outgoing spherical geometry computation. This was     #
! #       only correct for the upwelling case, and actually     #
! #       incorrect (in a small way) for the (rarely-used)      #
! #       downwelling case. Now, there are separate outgoing    #
! #       sphergeom routine for "_up" and "dn" and these are    #
! #       called separately in the various SSCORR_OUTGOING      #
! #       single scatter routines.                              #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_CHAPMAN                                 &
          ( DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,    & !  Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,         & !  Input
            EARTH_RADIUS, RFINDEX_PARAMETER,              & !  Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, & !  Input
            CHAPMAN_FACTORS, SZA_LOCAL_INPUT,             & !  Output
            FAIL, MESSAGE, TRACE )                          !  Output

!  The following options apply:

!   1. If the plane-parallel flag is on, no further inputs are required

!   2. If the plane-parallel flag is off, then Pseudo-spherical:
!       (a) Straight line geometry, must specify
!               Earth_radius, height grid
!       (b) Refractive geometry, must specify
!               Earth_radius, height grid
!               pressure grid, temperature grid

!  The logic will be checked before the module is called.

!  Newly programmed by R. Spurr, RT SOLUTIONS Inc. 5/5/05.

!    Based round a call to a pure geometry module which returns slant
!    path distances which was adapted for use in the Radiant model by
!    R. Spurr during an OCO L2 intensive April 24-29, 2005.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAXLAYERS

!  implicit none

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  flag for plane parallel case

      LOGICAL         :: DO_PLANE_PARALLEL

!  flag for refractive geometry

      LOGICAL         :: DO_REFRACTIVE_GEOMETRY

!  Number of layers

      INTEGER         :: NLAYERS

!  number of solar beams to be processed

      INTEGER         :: NBEAMS

!  number of fine layers within coarse layers

      INTEGER         :: FINEGRID(MAXLAYERS)

!  TOA solar zenith angles

      REAL(fpk)       :: BEAM_SZAS ( MAXBEAMS )

!  Earth radius (km)

      REAL(fpk)       :: EARTH_RADIUS
        
!  Refractive index parametaer (Born-Wolf approximation)

      REAL(fpk)       :: RFINDEX_PARAMETER

!  Coarse grids of heights, pressures and temperatures

      REAL(fpk)       :: HEIGHT_GRID     (0:MAXLAYERS)
      REAL(fpk)       :: PRESSURE_GRID   (0:MAXLAYERS)
      REAL(fpk)       :: TEMPERATURE_GRID(0:MAXLAYERS)

!  Output arguments
!  ----------------

!  Chapman factors

      REAL(fpk)       :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  solar zenith angles at nadir

      REAL(fpk)       :: SZA_LOCAL_INPUT(0:MAXLAYERS,MAXBEAMS)

!  output status

      LOGICAL         :: FAIL
      CHARACTER*(*)   :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  number of iterations (refractive case only)
!      This is debug output

      INTEGER         :: ITERSAVE(MAXLAYERS)

!  other local variables

      INTEGER      :: IB
      REAL(fpk)    :: SUN0

!  get spherical optical depths
!  ----------------------------

!  start beam loop

      DO IB = 1, NBEAMS

        SUN0 = BEAM_SZAS(IB)

        CALL BEAM_GEOMETRY_PREPARE                         &       
          ( MAXBEAMS, MAXLAYERS, IB, NLAYERS, FINEGRID,    & ! Input
            SUN0, EARTH_RADIUS, RFINDEX_PARAMETER,         & ! Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,     & ! Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,  & ! Input
            CHAPMAN_FACTORS, SZA_LOCAL_INPUT,             & ! Input
            ITERSAVE, FAIL, MESSAGE )                        ! Output

!  return if failed

        IF ( FAIL ) THEN
          TRACE = ' BEAM_GEOMETRY_PREPARE Geometry failure'
          RETURN
        ENDIF

!  end beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_CHAPMAN

!

SUBROUTINE BEAM_GEOMETRY_PREPARE                           & 
          ( MAXBEAMS, MAXLAYERS, IBEAM, NLAYERS, FINEGRID, & ! Input
            SZA_GEOM_TRUE, REARTH, RFINDEX_PARAMETER,      & ! Input
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,     & ! Input
            HEIGHTS, PRESSURES, TEMPERATURES,              & ! Input
            CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT,             & ! Input
            ITERSAVE, FAIL, MESSAGE )                        ! Output

!  Implicit none

      IMPLICIT NONE
        
!  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Coarse layering is input to the module. Values of Z, P, T  are
!  given at the layer boundaries, with the first value (index 0) at TOA.

!  The refractive geometry is assumed to start at the TOA level.

!  We also require the earth radius and the refractive index parameter
!   (For the Born-Wolf approximation)

!  There is no refraction if the flag DO_REFRACTIVE_GEOMETRY is not set.
!  In this case we do not require pressure and temperature information.
!  The calculation will then be for geometric rays.

!  The plane parallel Flag and the refractive geometry flag should not
!  both be true - this should be checked outside.

!  In the refracting case, fine-gridding of pressure and temperature is
!  done internally, temperature is interpolated linearly with height,
!  and pressure log-linearly. The refraction uses Snell's law rule.
!  Finelayer gridding assumes equidistant heights within coarse layers
!  but the number of fine layers can be varied

!  Output is specified at coarse layer boundaries

!  Module is stand-alone.

!  Reprogrammed for the OCO L2 algorithm
!   R. Spurr, RT Solutions, Inc.   April 27, 2005

!  Intended use in LIDORT and Radiant RT models.

!  Input arguments
!  ===============

!  input dimensioning

      INTEGER      :: MAXLAYERS, MAXBEAMS

!  Beam index

      INTEGER      :: IBEAM

!  number of coarse layers

      INTEGER      :: NLAYERS

!  number of fine layers within coarse layers

      INTEGER      :: FINEGRID(MAXLAYERS)

!  True solar zenith angle (degrees)

      REAL(fpk)    :: SZA_GEOM_TRUE

!  Earth radius (km)

      REAL(fpk)    :: REARTH
        
!  Refractive index parametaer (Born-Wolf approximation)

      REAL(fpk)    :: RFINDEX_PARAMETER

!  flag for plane parallel case

      LOGICAL      :: DO_PLANE_PARALLEL
        
!  flag for refractive geometry

      LOGICAL      :: DO_REFRACTIVE_GEOMETRY

!  Coarse grids of heights, pressures and temperatures

      REAL(fpk)    :: HEIGHTS     (0:MAXLAYERS)
      REAL(fpk)    :: PRESSURES   (0:MAXLAYERS)
      REAL(fpk)    :: TEMPERATURES(0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

      REAL(fpk)      :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  solar zenith angles at nadir

      REAL(fpk)      :: SZA_LEVEL_OUTPUT(0:MAXLAYERS,MAXBEAMS)

!  number of iterations (refractive case only)
!   This is debug output

      INTEGER         :: ITERSAVE(MAXLAYERS)

!  output status

      LOGICAL         :: FAIL
      CHARACTER*(*)   :: MESSAGE

!  Local variables
!  ===============

!  local dimensioning

      INTEGER, PARAMETER ::  LOCAL_MAXLAYERS     = 200
      INTEGER, PARAMETER ::  LOCAL_MAXFINELAYERS = 10
      
!  fine layer gridding for refraction

      REAL(fpk)    :: ZRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      REAL(fpk)    :: PRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)
      REAL(fpk)    :: TRFINE(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS)

!  local height arrays

      REAL(fpk)    :: H(0:LOCAL_MAXLAYERS)
      REAL(fpk)    :: DELZ(LOCAL_MAXLAYERS)

!  help variables

      INTEGER      :: N, J, NRFINE, K, ITER, MAXF,IB
      LOGICAL      :: LOOP
      REAL(fpk)    :: GM_TOA, TH_TOA, MU_TOA, MU_NEXT
      REAL(fpk)    :: Z1, Z0, Z, T1, T0, T, P1, P0, Q1, Q0, Q
      REAL(fpk)    :: FU, FL, DEG_TO_RAD

      REAL(fpk)    :: LAYER_DIST, MU_PREV, STH1, SINTH1, STH2, SINTH2, LOCAL_SUBTHICK, &
                      PHI, PHI_0, PHI_CUM, SINPHI, DELPHI, REFRAC, RATIO,              &
                      RE_LOWER, RE_UPPER, DIST, STH2D, SINTH2D, SNELL

!  Standard temperature (K) and pressure (mbar).

      REAL(fpk), PARAMETER  :: T_STANDARD = 273.16D0
      REAL(fpk), PARAMETER  :: P_STANDARD = 1013.25D0
      REAL(fpk), PARAMETER  :: STP_RATIO = T_STANDARD / P_STANDARD

!  Loschmidt's number (particles/cm2/km).

      REAL(fpk), PARAMETER  :: RHO_STANDARD = 2.68675D+24

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0,IB) = 0.0D0
      DO N = 1, NLAYERS
        SZA_LEVEL_OUTPUT(N,IB) = 0.0D0
        DO K = 1, NLAYERS
          CHAPMAN_FACTORS(N,K,IB) = 0.0D0
        ENDDO
      ENDDO

      FAIL    = .FALSE.
      MESSAGE = ' '

!  check local dimensioning

      IF ( LOCAL_MAXLAYERS .LT. NLAYERS ) THEN
        MESSAGE = 'local coarse layer dimensioning insufficient'
        FAIL = .TRUE.
        RETURN
      ENDIF

!  Check fine layers do not exceed local dimensions assigned

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
       MAXF = 0
       DO N = 1, NLAYERS
         MAXF = MAX(MAXF,FINEGRID(N))
       ENDDO
       IF ( LOCAL_MAXFINELAYERS .LT. MAXF ) THEN
         MESSAGE = 'local fine layer dimensioning insufficient'
         FAIL = .TRUE.
         RETURN
       ENDIF
      ENDIF
  
!  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

!  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )

!  initialize

      STH2D  = 0.0D0
        
!  derive the fine values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

        Z0 = HEIGHTS(0)
        P0 = PRESSURES(0)
        T0 = TEMPERATURES(0)
        Q0 = DLOG(P0)
        DO N = 1, NLAYERS
          NRFINE = FINEGRID(N)
          LOCAL_SUBTHICK = DELZ(N) / DBLE(NRFINE)
          P1 = PRESSURES(N)
          Z1 = HEIGHTS(N)
          T1 = TEMPERATURES(N)
          Q1 = DLOG(P1)
          ZRFINE(N,0) = Z0
          PRFINE(N,0) = P0
          TRFINE(N,0) = T0
          DO J = 1, NRFINE - 1
            Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
            FL = ( Z0 - Z ) / DELZ(N)
            FU = 1.0d0 - FL
            Q  = FL * Q1 + FU * Q0
            T  = FL * T0 + FU * T1
            PRFINE(N,J) = DEXP (  Q )
            TRFINE(N,J) = T
            ZRFINE(N,J) = Z
          ENDDO
          PRFINE(N,NRFINE) = P1
          TRFINE(N,NRFINE) = T1
          ZRFINE(N,NRFINE) = Z1
!              write(*,'(i3,11F10.4)')N,(PRFINE(N,J),J=0,NRFINE)
          Z0 = Z1
          P0 = P1
          T0 = T1
          Q0 = Q1
        ENDDO
      ENDIF

!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF
      
!  Refractive Geometry case
!  ========================

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN

!  Start value of SZA cosine

        MU_PREV = MU_TOA

!  start layer loop

        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1 = DASIN(SINTH1)
          PHI_0 = TH_TOA - STH1
          NRFINE = FINEGRID(N)

!  iteration loop

          ITER = 0
          LOOP = .TRUE.
          DO WHILE (LOOP.AND.ITER.LT.100)
            ITER = ITER + 1
            PHI_CUM = 0.0D0
            RE_UPPER = ZRFINE(1,0) + REARTH
            RATIO  = PRFINE(1,0) * STP_RATIO / TRFINE(1,0)
            REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
            SNELL = REFRAC * RE_UPPER * SINTH1
            DO K = 1, N
              LAYER_DIST = 0.0D0
              LOCAL_SUBTHICK = DELZ(K) / DBLE(NRFINE)
              DO J = 0, NRFINE - 1
                RATIO  = PRFINE(K,J) * STP_RATIO / TRFINE(K,J)
                REFRAC = 1.0D0 + RFINDEX_PARAMETER * RATIO
                RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                SINTH2 = SNELL/ (REFRAC * RE_UPPER )
                IF ( SINTH2.GT.1.0D0 ) SINTH2 = 1.0D0
                STH2 = DASIN(SINTH2)
                SINTH2D = RE_UPPER * SINTH2 / RE_LOWER
                IF ( SINTH2D .GT. 1.0D0 ) THEN
                  MESSAGE = 'refraction yields angles > 90 some levels'
                  FAIL = .TRUE.
                  RETURN
                ENDIF
                STH2D = DASIN(SINTH2D)
                PHI = STH2D - STH2
                SINPHI = DSIN(PHI)
                PHI_CUM = PHI_CUM + PHI
                DIST = RE_UPPER * SINPHI / SINTH2D
                LAYER_DIST = LAYER_DIST +  DIST
                RE_UPPER = RE_LOWER
              ENDDO
              CHAPMAN_FACTORS(N,K,IB) = LAYER_DIST / DELZ(K)
            ENDDO

!  examine convergence

            DELPHI = PHI_0 - PHI_CUM
            LOOP = (DABS(DELPHI/PHI_CUM).GT.1.0D-4)

!  Fudge factors to speed up the iteration

            IF ( SZA_GEOM_TRUE .GT. 88.7D0 ) THEN
              STH1 = STH1 + 0.1 * DELPHI
              PHI_0 = TH_TOA - STH1
            ELSE IF ( SZA_GEOM_TRUE .LT. 80.0D0 ) THEN
              PHI_0 = PHI_CUM
              STH1 = TH_TOA - PHI_0
            ELSE
              STH1 = STH1 + 0.3 * DELPHI
              PHI_0 = TH_TOA - STH1
            ENDIF
            SINTH1 = DSIN(STH1)

          ENDDO

!  failure

          IF ( LOOP ) THEN
            MESSAGE = 'refractive iteration not converged'
            FAIL = .TRUE.
            RETURN
          ENDIF

!  Update and save angle output

          MU_NEXT = DCOS(STH2D)
          MU_PREV = MU_NEXT
          SZA_LEVEL_OUTPUT(N,IB) = DACOS(MU_NEXT) / DEG_TO_RAD
          ITERSAVE(N) = ITER

        ENDDO

!  Straight line geometry
!  ======================

      ELSE
 
        DO N = 1, NLAYERS

!  start values

          SINTH1 = GM_TOA * H(N) / H(0)
          STH1   = DASIN(SINTH1)
          RE_UPPER = H(0)

!  solar zenith angles are all the same = input value

          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N

          DO K = 1, N

!  sine-rule; PHI = earth-centered angle

            RE_LOWER = RE_UPPER - DELZ(K)
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = DASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = DSIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

!  re-set

            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2

          ENDDO

!  finish main layer loop

        ENDDO

!  Finish

      ENDIF

!  end of routine

      RETURN
END SUBROUTINE BEAM_GEOMETRY_PREPARE

!

subroutine outgoing_sphergeom_fine_up                              &
       ( maxlayers, maxfine, maxpartials,                          & ! Input
         do_fine, do_partials, nlayers, nfine,                     & ! Input
         n_partials, partials_idx,                                 & ! Input
         heights, heights_p, eradius,                              & ! Input
         alpha_boa, theta_boa, phi_boa,                            & ! Input
         sunpaths,      radii,      ntraverse,      alpha_all,     & ! Output
         sunpaths_fine, ntraverse_fine, alpha_fine,                & ! Output
         sunpaths_p,    ntraverse_p,    alpha_p,                   & ! Output
         partials_fineidx, cosscat_up, fail, message )               ! Output

!  Completely stand-alone geometry routine for the outgoing correction
!     This is applicable to the Upwelling path geometry

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment

!  Version 3.3. September 2007, Partial layer geometries added
!  Version 3.5. Revision to separate upwelling/downwelling geometry routines

      implicit none

!  inputs

      integer       :: maxlayers, maxfine, maxpartials
      integer       :: nlayers, nfine
      logical       :: do_fine, do_partials
      integer       :: n_partials
      integer       :: partials_idx    (maxpartials)
      real(fpk)     :: eradius, heights (0:maxlayers)
      real(fpk)     :: heights_p(maxpartials)
      real(fpk)     :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      integer       :: ntraverse   (0:maxlayers)
      real(fpk)     :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk)     :: radii      (0:maxlayers)
      real(fpk)     :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      integer       :: ntraverse_fine(maxlayers,maxfine)
      real(fpk)     :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      real(fpk)     :: alpha_fine    (maxlayers,maxfine)

!  Partial layer output (geometry)

      integer       :: ntraverse_p(maxpartials)
      real(fpk)     :: sunpaths_p (maxpartials,maxlayers)
      real(fpk)     :: alpha_p    (maxpartials)

!  scattering angles

      real(fpk)     :: cosscat_up (0:maxlayers)

!   Local arrays

      real(fpk)     :: radii_fine    (maxlayers,maxfine)
      real(fpk)     :: radii_p    (maxpartials)
      real(fpk)     :: lospaths(maxlayers)

      real(fpk)     :: lospaths_p_up(maxpartials)
      integer       :: partials_fineidx(maxpartials)

      real(fpk)     :: theta_all  (0:maxlayers)
      real(fpk)     :: phi_all    (0:maxlayers)

!  Status output

      logical       :: fail
      character*(*) :: message

!  Local

      logical       :: direct_sun
      integer       :: n, j, k, krad, n1, ut, np
      real(fpk)     :: deg_to_rad, ex, ey, ez, px, py, pz
      real(fpk)     :: salpha_boa, calpha_boa, sphi_boa
      real(fpk)     :: stheta_boa, ctheta_boa, cphi_boa
      real(fpk)     :: ksi, cksi, sksi, xicum, tangr, fac
      real(fpk)     :: ctheta, stheta, calpha, salpha, cphi
      real(fpk)     :: b, sth0, th0, ks1, sth1, th1, theta_p

!  Local arrays associated with fine grid output

      integer, parameter :: maxlocalfine = 20
      logical            :: direct_sunf(maxlocalfine)
      real(fpk)          :: difz, dfine1, saf, xicum0, path
      real(fpk)          :: thetaf(maxlocalfine), xicumf, difa
      real(fpk)          :: cthetaf(maxlocalfine)
      real(fpk)          :: sthetaf(maxlocalfine)
      real(fpk)          :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90]); Check it!'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90]); Check it!'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif 

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo

!  Zero the partial paths, initialize the traverse number (nominal)

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          ntraverse_p(ut) = np
          do k = 1, np
            sunpaths_p(ut,k) = 0.0d0
          enddo
        enddo
      endif

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = 0.0d0
         enddo
        enddo
       enddo
      endif

!  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))

      cosscat_up (nlayers) = - calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif

!  Partial levels
!   Find the fine level just below partial point = (0,1,2,...nfine)

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          radii_p(ut) = eradius + heights_p(ut)
          j = 1
          do while (radii_p(ut).gt.radii_fine(np,j).and.j.le.nfine)
           j = j + 1
          enddo
          partials_fineidx(ut) = j - 1
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = 0.0d0
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = dsin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = dasin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = dsin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  PARTIAL LAYERS
!  --------------

        if ( do_partials ) then
          do ut = 1, n_partials
            np = partials_idx(ut)

!  Los angle and paths

            alpha_p(ut) = alpha_all(nlayers)
!            lospaths_p_dn(ut) = radii(np-1)-radii_p(ut)
            lospaths_p_up(ut) = radii_p(ut)-radii(np)

!  Overhead sun

            if (stheta_boa.eq.0.0d0 ) then
              do k = 1, np-1
                sunpaths_p(ut,k) = radii(k-1)-radii(k)
              enddo
              sunpaths_p(ut,np) = radii(np-1)-radii_p(ut)
            endif

!  Non-overhead sun

            if (stheta_boa.gt.0.0d0 ) then
              sth0 = stheta_boa
              th0  = theta_all(np)
              sth1 = sth0*radii_p(ut)/radii(np-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
              sth0 = sth1
              th0  = th1
              do k = np-1, 1, -1
                sth1 = sth0*radii(k)/radii(k-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
                sth0 = sth1
                th0  = th1
              enddo
            endif

!  Finish partial layer stuff

          enddo
        endif

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing spehricity geometry
!  ============================

!  define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.0.0d0 ) then
!        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * dsin(xicum)
!        py = 0.0d0
!        pz =   radii(0) * dcos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.0.0d0 ) then
!          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
!        endif
!      endif

!  initialise los cumulative angle

      xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1
  
!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = dsin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = dcos(thetaf(j))
             sthetaf(j) = dsqrt(1.0d0-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * dsin(xicumf)
             py = 0.0d0
             pz =   radii_fine(n1,j) * dcos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = dsqrt(1.0d0-cthetaf(j)*cthetaf(j))
             thetaf(j)  = dacos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = dsin(xicum)
!        f2y = 0.0d0
!        f2z = dcos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.1.0d0)  cphi = 1.0d0
!        if ( cphi.lt.-1.0d0) cphi = -1.0d0
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle

        cosscat_up(n) = cosscat_up(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!   Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,-1)
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)   

! DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1,-1)

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = dasin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

! DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = dasin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = dasin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = dasin(sth1)
            ks1 = th1-th0
            path = dsin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + 2.0d0*radii_fine(n1,j)*dcos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*dcos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tangent
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = dasin(sth1)
              ks1 = th1-th0
              path = path + 2.0d0 * dsin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1              
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = dasin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=2.0d0*radii(krad)*dcos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*dcos(th1)
            endif

! DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  PARTIALS Section
!  ----------------

!  Partial angles and lospaths

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)

!  establish the LOS angle, and lospaths

          th0 = alpha_all(np)
          sth0 = dsin(th0)
          salpha = radii(np)*sth0 / radii_p(ut)
          alpha_p(ut) = dasin(salpha)
          calpha = dsqrt(1.0d0-salpha*salpha)

          ksi = alpha_all(np) - alpha_p(ut)
          lospaths_p_up(ut) = dsin(ksi)*radii(np)/salpha

!          ksi = alpha_p(ut) - alpha_all(np-1)
!          lospaths_p_dn(ut) = dsin(ksi)*radii(np-1)/salpha

!  Check debug
!          lospaths_p_dn(ut) = lospaths(np)-lospaths_p_up(ut)
!          write(*,*)lospaths_p_up(ut),lospaths_p_dn(ut),
!     &              radii_p(ut)-radii(np)
!          pause

!  Establish the solar angle, both naidr =sun-at-BOA and generally

          xicum = alpha_all(nlayers) - alpha_p(ut)
          direct_sun = .true.
          if (stheta_boa.eq.0.0d0 ) then
            theta_p = xicum
            ctheta = dcos(theta_p)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
          else
            px = - radii_p(ut) * dsin(xicum)
            py = 0.0d0
            pz =   radii_p(ut) * dcos(xicum)
            b = ex*px + ey*py + ez*pz
            ctheta = -b/radii_p(ut)
            direct_sun = (direct_sun.and.ctheta.ge.0.d0)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            theta_p = dacos(ctheta)
          endif

!         write(*,*)direct_sun,theta_all(Np-1),theta_p,theta_all(np)

!  Sun paths, Direct sun at layer top
!   First calculate sunpath for layer containing partial point
!   Then work upwards to TOA using the sine rule.

          if ( direct_sun ) then
            sth0 = stheta
            th0  = theta_p
            sth1 = sth0*radii_p(ut)/radii(np-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
            sth0 = sth1
            th0  = th1
            do k = np-1, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
          endif

!  Debug
! 14       format(a11,1x,10f10.6)
!          j = partials_fineidx(ut)
!          write(*,*)ut,ntraverse_p(ut)
!          write(*,*)ut,j,np,lospaths_p_up(ut)/lospaths(np)
!           write(*,14)'regular feb',(sunpaths(np,k),k=np,1,-1)
!           do j = 1, nfine
!             write(*,14)'regular feb',(sunpaths_fine(np,k,j),k=np,1,-1)
!             if (j.eq.partials_fineidx(ut))
!     &        write(*,14)'regular ut ',(sunpaths_p(ut,k),k=np,1,-1)
!           enddo
!           write(*,14)'regular feb',0.0d0,(sunpaths(np-1,k),k=np-1,1,-1)
!          
!  Sun paths, Not direct sun , with tangent point. Code similar to above.
!  TANGR = tangent point radius for this ray
!   ntraverse_p(ut) = Number of layers traversed by ray.

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

          if (.not.direct_sun ) then
            tangr = stheta*radii(n)
            krad = nlayers
            do while (tangr.gt.radii(krad))
              krad = krad - 1
            enddo
            ntraverse_p(ut) = krad + 1
            sth0 = tangr/radii(0)
            th0 = dasin(sth0)
            do k = 1, krad
              sth1 = radii(k-1)*sth0/radii(k)
              th1 = dasin(sth1)
              ks1 = th1-th0
              fac = 1.0d0
              if ( k.gt.np) fac = 2.0d0
              sunpaths_p(ut,k) = fac*dsin(ks1)*radii(k-1)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            sunpaths_p(ut,krad+1)=2.0d0*radii_p(ut)*dcos(th1)   
          endif

!  Finish partials loop

        enddo
      endif

!  Finish

      return
end subroutine outgoing_sphergeom_fine_up

!

subroutine outgoing_sphergeom_fine_dn                              &
       ( maxlayers, maxfine, maxpartials,                          & ! Input
         do_fine, do_partials, nlayers, nfine,                     & ! Input
         n_partials, partials_idx,                                 & ! Input
         heights, heights_p, eradius,                              & ! Input
         alpha_boa, theta_boa, phi_boa,                            & ! Input
         sunpaths,      radii,      ntraverse,      alpha_all,     & ! Output
         sunpaths_fine, ntraverse_fine, alpha_fine,                & ! Output
         sunpaths_p,    ntraverse_p,    alpha_p,                   & ! Output
         partials_fineidx, cosscat_dn, fail, message )               ! Output

!  Completely stand-alone geometry routine for the outgoing correction
!    This is applicable to the Downwelling path geoemtry

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.3. September 2007, Partial layer geometries added

      implicit none

!  inputs

      integer       :: maxlayers, maxfine, maxpartials
      integer       :: nlayers, nfine
      logical       :: do_fine, do_partials
      integer       :: n_partials
      integer       :: partials_idx    (maxpartials)
      real(fpk)     :: eradius, heights (0:maxlayers)
      real(fpk)     :: heights_p(maxpartials)
      real(fpk)     :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      integer       :: ntraverse   (0:maxlayers)
      real(fpk)     :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk)     :: radii      (0:maxlayers)
      real(fpk)     :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      integer       :: ntraverse_fine(maxlayers,maxfine)
      real(fpk)     :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      real(fpk)     :: alpha_fine    (maxlayers,maxfine)

!  Partial layer output (geometry)

      integer       :: ntraverse_p(maxpartials)
      real(fpk)     :: sunpaths_p (maxpartials,maxlayers)
      real(fpk)     :: alpha_p    (maxpartials)

!  scattering angles

      real(fpk)     :: cosscat_dn (0:maxlayers)

!   Local arrays

      real(fpk)     :: radii_fine    (maxlayers,maxfine)
      real(fpk)     :: radii_p    (maxpartials)
      real(fpk)     :: lospaths(maxlayers)

      real(fpk)     :: lospaths_p_dn(maxpartials)
      integer       :: partials_fineidx(maxpartials)

      real(fpk)     :: theta_all  (0:maxlayers)
      real(fpk)     :: phi_all    (0:maxlayers)

!  Status output

      logical       :: fail
      character*(*) :: message

!  Local

      logical       :: direct_sun
      integer       :: n, j, k, krad, n1, ut, np
      real(fpk)     :: deg_to_rad, ex, ey, ez, px, py, pz
      real(fpk)     :: salpha_boa, calpha_boa, sphi_boa
      real(fpk)     :: stheta_boa, ctheta_boa, cphi_boa
      real(fpk)     :: ksi, cksi, sksi, xicum, tangr, fac
      real(fpk)     :: ctheta, stheta, calpha, salpha, cphi
      real(fpk)     :: b, sth0, th0, ks1, sth1, th1, theta_p

!  Local arrays associated with fine grid output

      integer, parameter :: maxlocalfine = 20
      logical            :: direct_sunf(maxlocalfine)
      real(fpk)          :: difz, dfine1, saf, xicum0, path
      real(fpk)          :: thetaf(maxlocalfine), xicumf, difa
      real(fpk)          :: cthetaf(maxlocalfine)
      real(fpk)          :: sthetaf(maxlocalfine)
      real(fpk)          :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90]); Check it!'
        fail    = .true.
        return
      endif
      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90]); Check it!'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif 

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo

!  Zero the partial paths, initialize the traverse number (nominal)

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          ntraverse_p(ut) = np
          do k = 1, np
            sunpaths_p(ut,k) = 0.0d0
          enddo
        enddo
      endif

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = 0.0d0
         enddo
        enddo
       enddo
      endif

!  start at BOA

      deg_to_rad = dacos(-1.0d0) / 180.0d0
      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))

      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa 

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif
 
!  Partial levels
!   Find the fine level just below partial point = (0,1,2,...nfine)
 
      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)
          radii_p(ut) = eradius + heights_p(ut)
          j = 1
          do while (radii_p(ut).gt.radii_fine(np,j).and.j.le.nfine)
           j = j + 1
          enddo
          partials_fineidx(ut) = j - 1
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = 0.0d0
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = dsin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = dasin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = dsin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  PARTIAL LAYERS
!  --------------

        if ( do_partials ) then
          do ut = 1, n_partials
            np = partials_idx(ut)

!  Los angle and paths

            alpha_p(ut) = alpha_all(nlayers)
            lospaths_p_dn(ut) = radii(np-1)-radii_p(ut)
!            lospaths_p_up(ut) = radii_p(ut)-radii(np)

!  Overhead sun

            if (stheta_boa.eq.0.0d0 ) then
              do k = 1, np-1
                sunpaths_p(ut,k) = radii(k-1)-radii(k)
              enddo
              sunpaths_p(ut,np) = radii(np-1)-radii_p(ut)
            endif

!  Non-overhead sun

            if (stheta_boa.gt.0.0d0 ) then
              sth0 = stheta_boa
              th0  = theta_all(np)
              sth1 = sth0*radii_p(ut)/radii(np-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
              sth0 = sth1
              th0  = th1
              do k = np-1, 1, -1
                sth1 = sth0*radii(k)/radii(k-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
                sth0 = sth1
                th0  = th1
              enddo
            endif

!  Finish partial layer stuff

          enddo
        endif

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing sphericity geometry
!  ============================

!  define Unit solar vector at BOA
!   Note change of sign for ex here.

      ex = + stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.0.0d0 ) then
!        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * dsin(xicum)
!        py = 0.0d0
!        pz =   radii(0) * dcos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.0.0d0 ) then
!          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
!        endif
!      endif

!  initialise los cumulative angle

      xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1
  
!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = dsin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = dcos(thetaf(j))
             sthetaf(j) = dsqrt(1.0d0-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * dsin(xicumf)
             py = 0.0d0
             pz =   radii_fine(n1,j) * dcos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = dsqrt(1.0d0-cthetaf(j)*cthetaf(j))
             thetaf(j)  = dacos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = dsin(xicum)
!        f2y = 0.0d0
!        f2z = dcos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.1.0d0)  cphi = 1.0d0
!        if ( cphi.lt.-1.0d0) cphi = -1.0d0
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle

        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_dn(n)-calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)
         phi_all(n)     = dacos(cphi)
        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!   Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,-1)
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)   

! DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1,-1)

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = dasin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

! DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = dasin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = dasin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = dasin(sth1)
            ks1 = th1-th0
            path = dsin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + 2.0d0*radii_fine(n1,j)*dcos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*dcos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tangent
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = dasin(sth1)
              ks1 = th1-th0
              path = path + 2.0d0 * dsin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1              
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = dasin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=2.0d0*radii(krad)*dcos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*dcos(th1)
            endif

! DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  PARTIALS Section
!  ----------------

!  Partial angles and lospaths

      if ( do_partials ) then
        do ut = 1, n_partials
          np = partials_idx(ut)

!  establish the LOS angle, and lospaths

          th0 = alpha_all(np)
          sth0 = dsin(th0)
          salpha = radii(np)*sth0 / radii_p(ut)
          alpha_p(ut) = dasin(salpha)
          calpha = dsqrt(1.0d0-salpha*salpha)

!          ksi = alpha_all(np) - alpha_p(ut)
!          lospaths_p_up(ut) = dsin(ksi)*radii(np)/salpha

          ksi = alpha_p(ut) - alpha_all(np-1)
          lospaths_p_dn(ut) = dsin(ksi)*radii(np-1)/salpha

!  Establish the solar angle, both naidr =sun-at-BOA and generally

          xicum = alpha_all(nlayers) - alpha_p(ut)
          direct_sun = .true.
          if (stheta_boa.eq.0.0d0 ) then
            theta_p = xicum
            ctheta = dcos(theta_p)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
          else
            px = - radii_p(ut) * dsin(xicum)
            py = 0.0d0
            pz =   radii_p(ut) * dcos(xicum)
            b = ex*px + ey*py + ez*pz
            ctheta = -b/radii_p(ut)
            direct_sun = (direct_sun.and.ctheta.ge.0.d0)
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            theta_p = dacos(ctheta)
          endif

!         write(*,*)direct_sun,theta_all(Np-1),theta_p,theta_all(np)

!  Sun paths, Direct sun at layer top
!   First calculate sunpath for layer containing partial point
!   Then work upwards to TOA using the sine rule.

          if ( direct_sun ) then
            sth0 = stheta
            th0  = theta_p
            sth1 = sth0*radii_p(ut)/radii(np-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_p(ut,np) = dsin(ks1)*radii_p(ut)/sth1
            sth0 = sth1
            th0  = th1
            do k = np-1, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths_p(ut,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
          endif

!  Debug
! 14       format(a11,1x,10f10.6)
!          j = partials_fineidx(ut)
!          write(*,*)ut,ntraverse_p(ut)
!          write(*,*)ut,j,np,lospaths_p_up(ut)/lospaths(np)
!           write(*,14)'regular feb',(sunpaths(np,k),k=np,1,-1)
!           do j = 1, nfine
!             write(*,14)'regular feb',(sunpaths_fine(np,k,j),k=np,1,-1)
!             if (j.eq.partials_fineidx(ut))
!     &        write(*,14)'regular ut ',(sunpaths_p(ut,k),k=np,1,-1)
!           enddo
!           write(*,14)'regular feb',0.0d0,(sunpaths(np-1,k),k=np-1,1,-1)
!          
!  Sun paths, Not direct sun , with tangent point. Code similar to above.
!  TANGR = tangent point radius for this ray
!   ntraverse_p(ut) = Number of layers traversed by ray.

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

          if (.not.direct_sun ) then
            tangr = stheta*radii(n)
            krad = nlayers
            do while (tangr.gt.radii(krad))
              krad = krad - 1
            enddo
            ntraverse_p(ut) = krad + 1
            sth0 = tangr/radii(0)
            th0 = dasin(sth0)
            do k = 1, krad
              sth1 = radii(k-1)*sth0/radii(k)
              th1 = dasin(sth1)
              ks1 = th1-th0
              fac = 1.0d0
              if ( k.gt.np) fac = 2.0d0
              sunpaths_p(ut,k) = fac*dsin(ks1)*radii(k-1)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            sunpaths_p(ut,krad+1)=2.0d0*radii_p(ut)*dcos(th1)   
          endif

!  Finish partials loop

        enddo
      endif

!  Finish

      return
end subroutine outgoing_sphergeom_fine_dn

!

subroutine multi_outgoing_adjustgeom                              &
         ( max_vza, max_sza, max_azm, n_vza, n_sza, n_azm,        & ! Input
          hsurface, eradius, do_adjust_surface,                   & ! Input
          alpha_boa, theta_boa, phi_boa,                          & ! Input
          alpha_ssa, theta_ssa, phi_ssa, fail, message, trace )     ! Output

! stand-alone geometry routine for adjusting the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height of new surface, earth radius.

!  Height grid here is artificial

      implicit none

!  inputs

      integer       :: max_vza, max_sza, max_azm
      integer       :: n_vza, n_sza, n_azm
      real(fpk)     :: eradius, hsurface
      logical       :: do_adjust_surface
      real(fpk)     :: alpha_boa (max_vza)
      real(fpk)     :: theta_boa (max_sza)
      real(fpk)     :: phi_boa   (max_azm)

!  outputs

      real(fpk)     :: alpha_ssa (max_vza)
      real(fpk)     :: theta_ssa (max_vza,max_sza,max_azm)
      real(fpk)     :: phi_ssa   (max_vza,max_sza,max_azm)
      logical       :: fail
      character*(*) :: message, trace

!  Local

      integer       :: j, i, k
      real(fpk)     :: deg_to_rad, ex, ey, ez, px, py, pz
      real(fpk)     :: salpha_boa, calpha_boa, sphi_boa
      real(fpk)     :: stheta_boa, ctheta_boa, cphi_boa
      real(fpk)     :: ksi, cksi, sksi, xicum, cosscat_up
      real(fpk)     :: phi_all, alpha_all, theta_all
      real(fpk)     :: ctheta, stheta, calpha, salpha, cphi
      real(fpk)     :: b,rssa
      character*2   :: c2

!  Initialise output

      fail = .false.
      message = ' '
      trace   = ' '

!  check range of inputs

      do j = 1, n_vza
       if ( alpha_boa(j).ge.90.0d0.or.alpha_boa(j).lt.0.0d0 ) then
        write(c2,'(I2)')J
        message = 'boa LOS angle outside range [0,90])'
        trace   = 'Change Boa Los angle, number '//C2
        fail    = .true.
        return
       endif
      enddo

      do k = 1, n_azm
        if ( phi_boa(k).lt.0.0d0   ) phi_boa(k) = - phi_boa(k)
        if ( phi_boa(k).gt.180.0d0 ) phi_boa(k) = 360.0d0 - phi_boa(k)
      enddo

      do i = 1, n_sza
       if ( theta_boa(i).ge.90.0d0.or.theta_boa(i).lt.0.0d0 ) then
        write(c2,'(I2)')I
        message = 'Boa SZA angle outside range [0,90])'
        trace   = 'Change Boa SZA angle, number '//C2
        fail    = .true.
        return
       endif
      enddo

!  No adjustment, just copy and exit

      if ( .not. do_adjust_surface ) then
        do j = 1, n_vza
          alpha_ssa(j)   = alpha_boa(j)
          do i = 1, n_sza
            do k = 1, n_azm
              theta_ssa(j,i,k) = theta_boa(i)
              phi_ssa(j,i,k)   = phi_boa(k)
            enddo
          enddo
        enddo
        return
      endif

!  conversion

      deg_to_rad = dacos(-1.0d0) / 180.0d0

!  Radius of surface

      rssa = hsurface + eradius

!  Start VZA loop

      do j = 1, n_vza

        alpha_all = alpha_boa(j) * deg_to_rad
        salpha_boa = dsin(alpha_all)
        calpha_boa = dcos(alpha_all)

!  Special case. Direct nadir viewing. Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

        if ( salpha_boa.eq.0.0d0 ) then
          alpha_ssa(j)   = alpha_boa(j)
          do i = 1, n_sza
            do k = 1, n_azm
              theta_ssa(j,i,k) = theta_boa(i)
              phi_ssa(j,i,k)   = phi_boa(k)
            enddo
          enddo
        endif

!  Los angle, general case

        if ( salpha_boa.ne.0.0d0 ) then

          salpha       = eradius * salpha_boa / rssa
          alpha_ssa(j) = dasin(salpha)
          calpha       = dcos(alpha_ssa(j))

!  Lospaths

          ksi  = alpha_all - alpha_ssa(j)
          sksi = dsin(ksi)
          cksi = dcos(ksi)
          xicum = ksi

!  output angle in degrees

          alpha_ssa(j) = alpha_ssa(j) / deg_to_rad

!  Los vector

          px = - rssa * dsin(xicum)
          py = 0.0d0
          pz =   rssa * dcos(xicum)

!  Start SZA loop

          do i = 1, n_sza

            theta_all = theta_boa(i) * deg_to_rad
            stheta_boa = dsin(theta_all)
            ctheta_boa = dcos(theta_all)

!  Start azimuth loop

            do k = 1, n_azm

              phi_all   = phi_boa(k)   * deg_to_rad
              cphi_boa  = dcos(phi_all)
              sphi_boa  = dsin(phi_all)

!  define Unit solar vector

              ex = - stheta_boa * cphi_boa
              ey = - stheta_boa * sphi_boa
              ez = - ctheta_boa
  
!  Sun angle

              b = ex*px + ey*py + ez*pz
              ctheta = -b/rssa
              stheta = dsqrt(1.0d0-ctheta*ctheta)
              theta_ssa(j,i,k) = dacos(ctheta)/deg_to_rad
              if ( ctheta.lt.0.0d0 ) then
                write(c2,'(I2)')J
                message = 'LOS-path SZA angle outside range [0,90])'
                trace   = 'Check inputs for LOS angle '//c2
                fail    = .true.
                return
              endif

!  scattering angle

              cosscat_up  = - calpha_boa * ctheta_boa + &
                              salpha_boa * stheta_boa * cphi_boa 

!  Fix phi by using constancy of scatter angle

              if ( phi_boa(k).eq.180.0d0 ) then
                phi_ssa(j,i,k) = phi_all
              else if ( phi_boa(k) .eq. 0.0d0 ) then
                phi_ssa(j,i,k) = 0.0d0
              else
                cphi = (cosscat_up+calpha*ctheta)/stheta/salpha
                if ( cphi.gt.1.0d0) cphi = 1.0d0
                if ( cphi.lt.-1.0d0) cphi = -1.0d0
                phi_ssa(j,i,k) = dacos(cphi)
              endif
              phi_ssa(j,i,k) = phi_ssa(j,i,k) / deg_to_rad

!  End Azimuth and solar loops

            enddo
          enddo

!  End of regular LOS case

        endif

!  End los loop

      enddo

!  Finish

      return
end subroutine multi_outgoing_adjustgeom

!

subroutine losonly_outgoing_adjustgeom &
       ( max_vza, n_vza, hsurface, eradius, do_adjust_surface, alpha_boa, &
         alpha_ssa, fail, message, trace )

      implicit none

!  Thermal only (no solar)

! Stand-alone geometry routine for adjusting the outgoing correction
!    starting inputs are the boa values of vza
!    need also the height of new surface, earth radius.

!  Height grid here is artificial

!  Inputs

      integer, intent(in) ::           max_vza
      integer, intent(in) ::           n_vza
      real(fpk), intent(in) ::         eradius, hsurface
      logical, intent(in) ::           do_adjust_surface
      real(fpk), intent(in) ::         alpha_boa (max_vza)

!  Outputs

      real(fpk), intent(out) ::        alpha_ssa (max_vza)
      logical, intent(out) ::          fail
      character (len=*), intent(inout) :: message, trace

!  Local

      integer ::    j
      real(fpk) ::  deg_to_rad, pie, pi2
      real(fpk) ::  salpha_boa, calpha_boa
      real(fpk) ::  alpha_all,  calpha, salpha, rssa
      character (len=2) :: c2

!  Initialise output

      fail    = .false.
      message = ' '
      trace   = ' '

!  Check range of inputs

      do j = 1, n_vza
       if ( alpha_boa(j).ge.90.0d0.or.alpha_boa(j).lt.0.0d0 ) then
        write(c2,'(i2)')j
        message = 'boa los angle outside range [0,90])'
        trace   = 'change boa los angle, number '//c2
        fail    = .true.
        return
       endif
      enddo

!  No adjustment, just copy and exit

      if ( .not. do_adjust_surface ) then
        do j = 1, n_vza
          alpha_ssa(j)   = alpha_boa(j)
        enddo
        return
      endif

!  Conversion

      pie = dacos(-1.0d0)
      pi2 = 2.0d0 * pie
      deg_to_rad = pie / 180.0d0

!  Radius of surface

      rssa = hsurface + eradius

!  Start vza loop

      do j = 1, n_vza

        alpha_all = alpha_boa(j) * deg_to_rad
        salpha_boa = dsin(alpha_all)
        calpha_boa = dcos(alpha_all)

!  Special case. direct nadir viewing. compute everything and exit.
!    (this is the same as the regular pseudo-spherical )

        if ( salpha_boa.eq.0.0d0 ) then
          alpha_ssa(j)   = alpha_boa(j)
        else
           salpha       = eradius * salpha_boa / rssa
           alpha_ssa(j) = dasin(salpha)
           calpha       = dcos(alpha_ssa(j))
           alpha_ssa(j) = alpha_ssa(j) / deg_to_rad
        endif

!  End los loop

      enddo

!  Finish

      return
end subroutine losonly_outgoing_adjustgeom

!  Finish

end module lidort_geometry

