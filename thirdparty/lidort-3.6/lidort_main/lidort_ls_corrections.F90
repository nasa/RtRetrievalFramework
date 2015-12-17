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
! #                                                         #
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

module lidort_ls_corrections

!  Parameter types

   USE LIDORT_PARS, only : fpk

public

! ###############################################################
! #                                                             #
! #      Version 3.3 outgoing Lambertian DB corrections         #
! #                                                             #
! #            LIDORT_LS_DBCORRECTION                           #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_LS_DBCORRECTION                                       &
        ( DO_SSCORR_OUTGOING, DO_UPWELLING, DO_BRDF_SURFACE,            & ! Input
          DO_REFLECTED_DIRECTBEAM, NLAYERS, NBEAMS, N_SURFACE_WFS,      & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP, FLUXMULT,            & ! Input
          LS_EXACTDB_BRDFUNC, ATTN_DB_SAVE,                             & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,      & ! Input
          SURFACEWF_DB )                                                  ! Output

!  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
!   Linearization with respect to LAMBERTIAN amplitude variable only !

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXLAYERS, MAX_PARTLAYERS, MAX_GEOMETRIES, MAXBEAMS, &
                              MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_USER_RELAZMS, &
                              MAX_SURFACEWFS, ZERO

      IMPLICIT NONE

!  Input Arguments
!  ---------------

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING

!  Flag for outgoing single scatter correction

      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING

!  BRDF flag

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Direct beam reflectance

      LOGICAL  , intent(in)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  Number of surface weighting functions

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Numbers

      INTEGER  , intent(in)  :: N_USER_RELAZMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Offsets for geometry indexing

      INTEGER  , intent(in)  :: N_GEOMETRIES
      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

! FLux

      REAL(fpk), intent(in)  :: FLUXMULT

!  Linearized Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), intent(in)  :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Atmospheric attenuation before reflection

      REAL(fpk), intent(in)  ::  ATTN_DB_SAVE(MAX_GEOMETRIES)

!  output optical depth masks and indices
!  off-grid optical depths (values, masks, indices)

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Outgoing sphericity Whole layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  Outgoing sphericity Partial-layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Output arguments
!  ----------------

!  Direct Bounce Surface linearization

      REAL(fpk), intent(inout) :: SURFACEWF_DB (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Local variables
!  ---------------

!  Linearized Cumulative Exact direct beam source terms

      REAL(fpk)  :: LS_DB_CUMSOURCE(MAX_SURFACEWFS, MAX_GEOMETRIES)

!  Help

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UT, UTA, UM, UA, NC, IB, V, K
      REAL(fpk)  :: TR, BRDF

!  first stage
!  -----------

!  Initialize

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_USER_LEVELS
          DO K = 1, N_SURFACE_WFS
            SURFACEWF_DB(K,UTA,V)  = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

!  initialize cumulative linearized source term, LAMBERTIAN
!    In this case, just original source term / albedo
!      ( as dependence on Lambertin albedo is linear )

      IF ( .not. DO_BRDF_SURFACE ) THEN
        DO V = 1, N_GEOMETRIES
          LS_DB_CUMSOURCE(1,V) = FLUXMULT * ATTN_DB_SAVE(V)
        ENDDO
      ENDIF

!  initialize cumulative linearized source term, BRDF

      IF ( DO_BRDF_SURFACE ) THEN
        DO K = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            DO IB = 1, NBEAMS
              IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  BRDF = FLUXMULT * LS_EXACTDB_BRDFUNC(K,UM,UA,IB)
                  LS_DB_CUMSOURCE(K,V) = BRDF * ATTN_DB_SAVE(V)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Transmittance of source term: upwelling recursion
!  -------------------------------------------------

!  initialize optical depth loop

      NC =  0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    loop over layers working upwards to level NUT

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS(N,V)
                ELSE
                  TR  = T_DELT_USERM(N,UM)
                ENDIF
                DO K = 1, N_SURFACE_WFS
                  LS_DB_CUMSOURCE(K,V) = TR * LS_DB_CUMSOURCE(K,V)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Offgrid output
!  -------------- 

!  Require partial layer transmittance

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR  = UP_LOSTRANS_UT(UT,V)
                ELSE
                  TR  = T_UTUP_USERM(UT,UM)
                ENDIF
                DO K = 1, N_SURFACE_WFS
                  SURFACEWF_DB(K,UTA,V) = TR * LS_DB_CUMSOURCE(K,V)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Ongrid output : Set final cumulative source directly

        ELSE
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                DO K = 1, N_SURFACE_WFS
                  SURFACEWF_DB(K,UTA,V) = LS_DB_CUMSOURCE(K,V)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Check for updating the recursion 

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end optical depth loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_LS_DBCORRECTION

!  end

end module lidort_ls_corrections

