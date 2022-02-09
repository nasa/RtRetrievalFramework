! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ############################################################### ! Old
! # Subroutines in this (former) LIDORT_corrections Module      # ! Old
! #            LIDORT_SSCORR_NADIR (master)                     # ! Old
! #      Version 3.2 outgoing sphericity correction             # ! Old
! #              3.3.  partial-layer integration                # ! Old
! #            LIDORT_SSCORR_OUTGOING (master)                  # ! Old
! #                 OUTGOING_INTEGRATION_UP                     # ! Old
! #                 OUTGOING_INTEGRATION_DN                     # ! Old
! #      Version 3.3 and beyond,DB correction                   # ! Old
! #            LIDORT_DBCORRECTION                              # ! Old
! ############################################################### ! Old

! ###############################################################
! #                                                             #
! #            SFO_MASTER_INTERFACE                             #
! #                                                             #
! ###############################################################

!  this is the standard interface to the FO code - NEW for LIDORT 3.7 !!!!

!  Version 3.0 - 3.7. Internal SSCORR/DBCORR routines (Old)
!  Version 3.8. Internal SSCORR/DBCORR routines removed.
!  Version 3.8.1, upgrade, June 2019
!    -- 4/9/19. Add FO Surface-leaving assignation, saved cumulative transmittance

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_UP, THETA_ALL, ALPHA (upwelling)
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_DN, THETA_ALL, ALPHA (downwelling)

!  2/28/21. Version 3.8.3. Add DO_DOUBLET Geometry flag, plus offset input arguments

      MODULE lidort_sfo_interface_m

      USE SFO_Master_m

      PRIVATE
      PUBLIC :: SFO_MASTER_INTERFACE

      CONTAINS

!

      SUBROUTINE SFO_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_PHASFUNC,                & ! Input SS control flags
        DO_DELTAM, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS, DO_OBSGEOM, DO_DOUBLET,       & ! Input RT Control flags
        DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,             & ! Input Optical/Surface flags
        NLAYERS, NFINELAYERS, NLEGEN_MOMENTS_INPUT, ND_OFFSET, NV_OFFSET, NA_OFFSET,        & ! Input numbers
        NBEAMS, BEAM_SZAS, N_USER_STREAMS, USER_STREAMS, N_USER_RELAZMS, USER_RELAZMS,      & ! Input geometry
        N_USER_LEVELS, levelmask_UP, levelmask_DN, N_PARTLAYERS,                            & ! Input levels  control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
        EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER,                                      & ! Input Flux/Heights
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
        DELTAU_VERT, PHASFUNC_UP, PHASFUNC_DN, TRUNC_FACTOR, BB_INPUT,                      & ! Inputs (Optical - Regular)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                        & ! Inputs (Optical - Surface)
        SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
        FO_INTENSITY_SS, FO_INTENSITY_DB, FO_INTENSITY_DTA, FO_INTENSITY_DTS,               & ! Output Intensity
        FO_INTENSITY_ATMOS, FO_INTENSITY_SURF, FO_INTENSITY,                                & ! Output Intensity
        CUMTRANS, LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA, SLTERM,                       & ! Output Auxiliary
        FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Exception handling

!  This routine is completely new for Version 3.8 of LIDORT, 3/3/17
!    Modeled after the similar routine in VLIDORT 2.8

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  2/28/21. Version 3.8.3. Add Doublet geometry flag and include offsets

      USE LIDORT_PARS_m, Only : MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS,          &
                                MAXLAYERS, MAXMOMENTS_INPUT, MAXMOMENTS, MAX_GEOMETRIES, MAXFINELAYERS, &
                                MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO, PIE, DEG_TO_RAD, fpk

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags. Water-leaving flag introduced for Version 3.8.1. 4/9/19.

      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) :: DO_SURFACE_EMISSION

      LOGICAL, INTENT(IN) :: DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) :: DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) :: DO_SSCORR_OUTGOING

      LOGICAL, INTENT(IN) :: DO_DELTAM
      LOGICAL, INTENT(IN) :: DO_UPWELLING
      LOGICAL, INTENT(IN) :: DO_DNWELLING
      LOGICAL, INTENT(IN) :: DO_PHASFUNC

      LOGICAL, INTENT(IN) :: DO_BRDF_SURFACE

!  2/28/21. Version 3.8.3. Add Doublet Geometry option

      LOGICAL, INTENT(IN) :: DO_OBSGEOM
      LOGICAL, INTENT(IN) :: DO_DOUBLET

      LOGICAL, INTENT(IN) :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) :: DO_WATER_LEAVING
      LOGICAL, INTENT(IN) :: DO_SL_ISOTROPIC

!  Control integers

      INTEGER, INTENT(IN) :: NLAYERS
      INTEGER, INTENT(IN) :: NFINELAYERS
      INTEGER, INTENT(IN) :: NLEGEN_MOMENTS_INPUT

!  Geometry

      INTEGER  , INTENT(IN) :: NBEAMS
      Real(fpk), INTENT(IN) :: BEAM_SZAS ( MAXBEAMS )
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      Real(fpk), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      INTEGER  , INTENT(IN) :: N_USER_RELAZMS
      Real(fpk), INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )

!  2/28/21. Version 3.8.3. Add OFFSETS inputs (Doublet ND and lattice NV/NA)

      INTEGER, INTENT (IN) :: ND_OFFSET ( MAXBEAMS )
      INTEGER, INTENT (IN) :: NV_OFFSET ( MAXBEAMS )
      INTEGER, INTENT (IN) :: NA_OFFSET ( MAXBEAMS, MAX_USER_STREAMS )

!  require the Level Mask inputs

      INTEGER, INTENT (IN) :: N_USER_LEVELS
      INTEGER, INTENT (IN) :: levelmask_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: levelmask_DN  ( MAX_USER_LEVELS )

!      Real(fpk), INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs.

      LOGICAL  , INTENT (IN) :: DO_PARTLAYERS
      INTEGER  , INTENT (IN) :: N_PARTLAYERS
      LOGICAL  , INTENT (IN) :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER  , INTENT (IN) :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER  , INTENT (IN) :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      Real(fpk), INTENT (IN) :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  other inputs

      Real(fpk), INTENT(IN) :: EARTH_RADIUS
      Real(fpk), INTENT(IN) :: HEIGHT_GRID ( 0:MAXLAYERS )

      Real(fpk), INTENT(IN) :: SS_FLUX_MULTIPLIER

!  Optical unscaled

      Real(fpk), INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )
      Real(fpk), INTENT(IN) :: BB_INPUT ( 0:MAXLAYERS )

!  optical scaled

      Real(fpk), INTENT(IN) :: DELTAU_VERT  ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: TRUNC_FACTOR ( MAXLAYERS )

!  Phase functions, introduced for Version 3.8 of LIDORT, and Version 1.5 of FO

      Real(fpk), INTENT(IN) :: PHASFUNC_UP ( MAX_GEOMETRIES, MAXLAYERS )
      Real(fpk), INTENT(IN) :: PHASFUNC_DN ( MAX_GEOMETRIES, MAXLAYERS )

!  surface BRDF and EMiss.

      Real(fpk), INTENT(IN) :: LAMBERTIAN_ALBEDO
      Real(fpk), INTENT(IN) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      Real(fpk), INTENT(IN) :: SURFBB
      Real(fpk), INTENT(IN) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Surface leaving inputs, 8/3/16, for Version 2.8

      Real(fpk), INTENT(IN) :: SLTERM_ISOTROPIC  ( MAXBEAMS )
      Real(fpk), INTENT(IN) :: SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Outputs
!  =======

!  Solar

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Thermal

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Composite

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  4/9/19. Additional output for the sleave correction

      real(fpk), Intent(out)     :: CUMTRANS ( max_user_levels, max_geometries )

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list

      real(fpk), Intent(out)     :: LOSTRANS_UP ( max_geometries, maxlayers )
      real(fpk), Intent(out)     :: LOSTRANS_DN ( max_geometries, maxlayers )

      real(fpk), Intent(out)     :: theta_all ( 0:maxlayers, max_geometries )
      real(fpk), Intent(out)     :: alpha     ( 0:maxlayers, max_geometries )

!  4/9/19. Surface leaving FO assignation

      Real(fpk), intent(out)     :: SLTERM ( MAX_GEOMETRIES )

!  Exception handling

      LOGICAL, INTENT (OUT)             :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1, TRACE_2

!  Local variables
!  ===============

!  Max dimensions
!  --------------

      INTEGER   :: MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, MAXFINE

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

      INTEGER   :: NGEOMS, NSZAS, NVZAS, NAZMS, NFINE, NMOMENTS_INPUT

!  Configuration inputs
!  --------------------

!  Flags (sphericity flags should be mutually exclusive)

      LOGICAL   :: DO_PLANPAR
      LOGICAL   :: DO_ENHANCED_PS

!  Lambertian surface flag

      LOGICAL   :: DO_LAMBERTIAN

!  General inputs
!  --------------

!  DTR = degrees-to-Radians

      Real(fpk) :: DTR

!  Critical adjustment for cloud layers

      LOGICAL   :: DoCrit
      Real(fpk) :: Acrit

!  Earth radius + heights. Partial heights.

      Real(fpk) :: ERADIUS
      Real(fpk) :: HEIGHTS ( 0:MAXLAYERS )
      Real(fpk) :: PARTIAL_HEIGHTS ( MAX_PARTLAYERS )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees). Note the dimensioning now,

      Real(fpk) :: obsgeom_boa  ( MAX_GEOMETRIES, 3 )
      Real(fpk) :: theta_boa    ( MAXBEAMS )             !SZA
      Real(fpk) :: alpha_boa    ( MAX_USER_STREAMS  )    !UZA
      Real(fpk) :: phi_boa      ( MAX_USER_RELAZMS  )    !RAA

!  Optical inputs
!  --------------

!  Solar flux

      Real(fpk) :: FLUX

!  Atmosphere

      Real(fpk) :: EXTINCTION  ( MAXLAYERS )
      Real(fpk) :: DELTAUS     ( MAXLAYERS )
      Real(fpk) :: OMEGA       ( MAXLAYERS )
      Real(fpk) :: PHASMOMS    ( MAXLAYERS, 0:MAXMOMENTS_INPUT )

!  For TMS correction

      Real(fpk) :: TRUNCFAC    ( MAXLAYERS )

!  Thermal inputs. Use directly, don't need to copy
!      Real(fpk) :: BB_INPUT ( 0:MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      Real(fpk) :: REFLEC ( MAX_GEOMETRIES )

!  Surface properties - emissive

      Real(fpk) :: EMISS ( MAX_USER_STREAMS )

!  Other variables
!  Other variables
!    -- 2/28/21. Version 3.8.3. Offsets now input

      integer   :: ns, nv, na
      INTEGER   :: G, L, N, UT
      INTEGER   :: LUM=1, LUA=1

!  Define SFO_MASTER inputs
!  ========================

!  Note: argument passing for variables defined elsewhere commented out here

!  Max dimensions.

      maxgeoms          = MAX_GEOMETRIES
      maxfine           = MAXFINELAYERS
      maxszas           = MAXBEAMS
      maxvzas           = MAX_USER_STREAMS
      maxazms           = MAX_USER_RELAZMS

!  Proxies
!     -- 2/28/21. Version 3.8.3. Add Double toption

      if ( do_ObsGeom ) then
         ngeoms  = NBEAMS
      else if ( do_Doublet ) then
         ngeoms  = NBEAMS * N_USER_STREAMS
      else
         ngeoms  = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
      endif
      nszas   = NBEAMS
      nvzas   = N_USER_STREAMS
      nazms   = N_USER_RELAZMS
      nfine           = NFINELAYERS
      nmoments_input  = NLEGEN_MOMENTS_INPUT

!  Offsets. -- 2/28/21. Version 3.8.3. Removed from here. Now inputs

!      na_offset = 0 ; nv_offset = 0
!      if ( .not. do_obsgeom ) then
!        do ns = 1, nszas
!          nv_offset(ns) = nvzas * nazms * (ns - 1) 
!          do nv = 1, nvzas
!            na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
!          enddo
!        enddo
!      endif

!  Configuration inputs. Removed "regular_ps" flag 9/17/16 (Redundant)

      if (DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (DO_SSCORR_NADIR) then
          do_enhanced_ps = .FALSE.
        else if (DO_SSCORR_OUTGOING) then
          do_enhanced_ps = .TRUE.
        end if
      end if
      do_lambertian      = .not. DO_BRDF_SURFACE

!  General inputs

      dtr    = DEG_TO_RAD

!      DoCrit = .FALSE. !set for now
!      Acrit  = ZERO    !set for now

!  Rob Fix 3/16/15. Default values changed 
!   [ Note for future - these should be inputs ]

      DoCrit = .TRUE.     
      Acrit  = 1.0d-12

!  Earth Radius and heights. Partials added, 9/17/16

      eradius             = EARTH_RADIUS
      heights(0:nlayers)  = HEIGHT_GRID(0:NLAYERS)
      do ut = 1, n_partlayers
        n = partlayers_layeridx(ut)
        partial_heights(ut) = heights(n-1) - partlayers_values(ut) * ( heights(n-1) - heights(n) )
      enddo

!  This code no longer required, now using same arrays as VLIDORT
!      fo_user_levels(1:n_user_levels) = INT(user_levels(1:n_user_levels))

!  Geometry inputs. all angles in Degrees.
!  These are good for Lattice/Doublet and ObsGeom (2/28/21. Version 3.8.3)

      obsgeom_boa = zero
      if ( do_ObsGeom ) then
         obsgeom_boa(1:ngeoms,1) = BEAM_SZAS     (1:ngeoms)
         obsgeom_boa(1:ngeoms,2) = USER_STREAMS (1:ngeoms)
         obsgeom_boa(1:ngeoms,3) = USER_RELAZMS  (1:ngeoms)
      endif
      theta_boa(1:NBEAMS)         = BEAM_SZAS    (1:NBEAMS)
      alpha_boa(1:N_USER_STREAMS) = USER_STREAMS (1:N_USER_STREAMS)
      phi_boa  (1:N_USER_RELAZMS) = USER_RELAZMS (1:N_USER_RELAZMS)

!  Optical inputs

      flux     = SS_FLUX_MULTIPLIER !=FLUX_FACTOR/(4.0d0*PIE)

! @@ Rob 7/30/13
!   If deltam-scaling, must use scaled deltau !!!

      extinction = 0.0d0
      if (DO_DELTAM) then
        extinction(1:nlayers) = DELTAU_VERT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT(1:NLAYERS)
      else
        extinction(1:nlayers) = DELTAU_VERT_INPUT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT_INPUT(1:NLAYERS)
      endif

!  Single scattering albedo is unscaled

      omega(1:nlayers)      = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Version 3.8 (FO 1.5). Optional use of Phase function input
!    Phasmoms will not be needed in this case

      Phasmoms = zero 
      if ( .not. do_phasfunc ) then
         do l=0,nmoments_input
            do n=1,nlayers
               Phasmoms(n,l) = PHASMOMS_TOTAL_INPUT(L,N)
            end do
         end do
      endif

! For TMS correction only

      truncfac = zero
      if (DO_DELTAM) then
        do n=1,nlayers
          truncfac(n) = TRUNC_FACTOR(N)
        end do
      end if

!  copy thermal BB input. Now used directly.
!      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

!  BRDF reflection now also for Lattice.
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option for BRDF reflection.

      reflec = ZERO
      if (do_lambertian) then
        do g = 1,ngeoms
          reflec(g) = LAMBERTIAN_ALBEDO
        end do
      else
        if ( do_ObsGeom ) then
           reflec(1:ngeoms) = EXACTDB_BRDFUNC(lum,lua,1:ngeoms)
        else  if ( do_doublet ) then
          do ns = 1, nszas ; do nv = 1, nvzas
            g = nd_offset(ns) + nv
            reflec(g) = EXACTDB_BRDFUNC(nv,LUA,ns)
          end do ; end do
        else
          do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
            g = na_offset(ns,nv) + na
            reflec(g) = EXACTDB_BRDFUNC(nv,na,ns)
          enddo ; enddo ; enddo
        end if
      end if

!  Surface leaving. 4/9/19 ==> This is now an output.
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 2.8.2). Add Do_Doublet option for surface leaving.

      SLTERM = zero
      if ( do_surface_leaving ) then
        if ( do_ObsGeom ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            SLTERM(1:ngeoms) = SLTERM_ISOTROPIC(1:ngeoms)
          ELSE
            SLTERM(1:ngeoms) = SLTERM_USERANGLES(LUM,LUA,1:ngeoms)
          ENDIF
        else if ( do_Doublet ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(g) = SLTERM_ISOTROPIC(ns)
            end do ; end do 
          ELSE
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(g) = SLTERM_USERANGLES(nv,LUA,ns)
            end do ; end do
          ENDIF
        else
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(g) = SLTERM_ISOTROPIC(ns)
            end do ; end do ; end do 
         ELSE
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(g) = SLTERM_USERANGLES(nv,na,ns)
            end do ; end do ; end do 
          ENDIF
        endif
      endif

!  Emissivity.

      emiss = zero
      if ( do_surface_emission ) then
        if ( do_ObsGeom ) then
          do g=1,ngeoms
            emiss(g) = USER_EMISSIVITY(g)
          end do
        else
          do nv=1,nvzas
            emiss(nv) = USER_EMISSIVITY(nv)
          end do
        end if
      endif

!  Call SFO_MASTER
!  ===============

!  Upgraded to FO Version 1.5, 7/11/16 in VLIDORT, 3/3/17 here in LIDORT
!   4/9/19. CUMTRANS now an additional output. add water-leaving control.
      
!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  2/28/21. Version 3.8.3. Other changes
!     -- Add Doublet flag to input list, add geometry offsets as inputs

      CALL SFO_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms,                                             & ! Input max dims
         maxlayers, max_partlayers, maxfine, maxmoments_input, max_user_levels,           & ! Input max dims
         do_solar_sources, do_thermal_emission, do_surface_emission, do_phasfunc,         & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_deltam, do_obsgeom, do_doublet,                   & ! Input flags (general)
         do_surface_leaving, do_water_leaving, do_partlayers, do_planpar, do_enhanced_ps, & ! Input flags (surface/geoms)
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,                     & ! Inputs (control-numbers)
         nd_offset, nv_offset, na_offset, n_user_levels, levelmask_up, levelmask_dn,      & ! Inputs (control-levels)
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx,      & ! Inputs (control-partial)
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                      & ! Input general
         obsgeom_boa, alpha_boa, theta_boa, phi_boa, flux,                                & ! Input geometry/Flux
         extinction, deltaus, omega, phasmoms, phasfunc_up, phasfunc_dn,                  & ! Input atmos optical
         truncfac, bb_input, surfbb, emiss, reflec, slterm,                               & ! Input thermal/surf optical
         fo_intensity_ss, fo_intensity_db, fo_intensity_dta, fo_intensity_dts,            & ! Output
         fo_intensity_atmos, fo_intensity_surf, fo_intensity,                             & ! Output
         cumtrans, lostrans_up, lostrans_dn, theta_all, alpha,                            & ! Auxiliary Output
         fail, message, trace_1, trace_2 )                                                  ! Exception-Handling

!  Done

      RETURN
END SUBROUTINE SFO_MASTER_INTERFACE

!  finish module

END MODULE lidort_sfo_interface_m


