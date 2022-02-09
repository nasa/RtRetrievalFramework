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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #  For linearizations involving atmospheric parameters        #
! #            LP_BVP_SOLUTION_MASTER                           #
! #            LP_BVP_COLUMN_SETUP                              #
! #            LP_BVP_SURFACE_SETUP                             #
! #            LP_BEAMSOLUTION_NEQK                             #
! #            LP_BEAMSOLUTION_NNEK                             #
! #                                                             #
! #  For linearizations (atmospheric) using telescoped BVP      #
! #            LP_BVPTEL_SOLUTION_MASTER                        #
! #            LP_BVPTEL_COLUMN_SETUP                           #
! #            LP_BVPTEL_SURFACE_SETUP                          #
! #                                                             #
! ###############################################################

!  2/28/21. Version 3.8.3.  Two Changes
!    -- BRDF Fourier inputs are defined locally for each Fourier component
!    -- Add CFUNC/DFUNC inpouts to LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER.

module lidort_lp_bvproblem_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

!  LIDORT Use dependencies

   USE lidort_aux_m,    only : DGBTRS, DGETRS
   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

private :: LP_BVP_COLUMN_SETUP, LP_BVP_SURFACE_SETUP, LP_BEAMSOLUTION_NEQK, LP_BEAMSOLUTION_NNEK, &
           LP_BVPTEL_COLUMN_SETUP, LP_BVPTEL_SURFACE_SETUP
public  :: LP_BVP_SOLUTION_MASTER, LP_BVPTEL_SOLUTION_MASTER

contains

SUBROUTINE LP_BVP_SOLUTION_MASTER &
      ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
        DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input
        DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, TAYLOR_ORDER,         & ! Input, Flags and order
        FOURIER, IBEAM, NLAYERS, NSTREAMS, NSTREAMS_2, NTOTAL,        & ! Input
        N_SUBDIAG, N_SUPDIAG, VARIATION_INDEX, N_WEIGHTFUNCS,         & ! Input
        DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
        SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
        DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,         & ! Input, BVP Bandmat
        T_DELT_EIGEN, XPOS, XNEG, LCON, MCON, ATERM_SAVE, BTERM_SAVE, & ! Input, Homogeneous/Greens
        CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,           & ! Input, Greens Function
        LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
        LP_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homogeneous solution 
        L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + Thermal
        NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,         & ! Output - Linearized Constants + PI
        STATUS, MESSAGE, TRACE )                                        ! Output - Exception handling

!  Linearization of the Boundary Problem Solution

!  2/28/21. Version 3.8.3.  BRDF Fourier inputs are defined locally for each Fourier component
!    -- Add CFUNC/DFUNC to Inputs. Rearrange I/O list slightly.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE
      LOGICAL  , intent(in)  ::  DO_WATER_LEAVING

!  Flag

      LOGICAL  , intent(in)  ::  DO_PLANE_PARALLEL

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)   :: TAYLOR_ORDER

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of streams

      INTEGER  , intent(in)   :: NSTREAMS, NSTREAMS_2

!  Number of layers

      INTEGER  , intent(in)   :: NLAYERS

!  BVProblem Band matrix control

      INTEGER  , intent(in)   :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Linearization control

      INTEGER  , intent(in)   :: VARIATION_INDEX, N_WEIGHTFUNCS

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  surface stuff
!  -------------

!  surface factor = 1+delta(m,0). Albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF_F defined locally (remove MAXMOMENTS dimension)

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )


!  4/9/19  RF_DIRECT_BEAM is the reflected beam (excludes SLTERM)

      REAL(fpk), intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SLTERM         ( MAXSTREAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  ::  IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  ::  SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(in)  ::  LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  ::  MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Homogeneous and thermal solution variables
!  ------------------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

      REAL(fpk), intent(in)  :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaing contribution

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_TRANS_ATMOS    ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  boundary condition flags

      LOGICAL    :: MODIFIED_BCL3, MODIFIED_BCL4

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI, CN

!  Other local help variables 

      INTEGER     :: I, Q, N, AA
      INTEGER     :: NS1, NS2         !mick eff 3/22/2017

!  Initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!mick eff 3/22/2017
!  Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Linearization of the regular BVP case
!  =====================================

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( VARIATION_INDEX .EQ. 1 )
      MODIFIED_BCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

!  Column vector setup
!    -- 2/28/21. Version 3.8.3. BRDF_F defined locally (remove MAXMOMENTS dimension)
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC arguments to input list

      CALL LP_BVP_COLUMN_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
           DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input, Flags
           DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, TAYLOR_ORDER,         & ! Input, Flag and order
           FOURIER, IBEAM, NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL,        & ! Input, Numbers
           MODIFIED_BCL3, MODIFIED_BCL4, VARIATION_INDEX, N_WEIGHTFUNCS, & ! Input
           DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
           SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
           DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
           T_DELT_EIGEN, XPOS, XNEG, LCON, MCON, ATERM_SAVE, BTERM_SAVE, & ! Input, Homogeneous/Greens
           CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,           & ! Input, Greens Function
           LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
           LP_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homogeneous
           L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
           L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_WEIGHTFUNCS, &
             BANDMAT2, MAXBANDTOTAL, IPIVOT,  COL2_WF, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) VARIATION_INDEX
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for layer '//CN//'DGBTRS call in LP_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, N_WEIGHTFUNCS
              NCON(I,N,Q) = COL2_WF(LCONMASK(I,N),Q)
              PCON(I,N,Q) = COL2_WF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
       ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, N_WEIGHTFUNCS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_WF, MAXSTREAMS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'Atmos_Wfs for 1-layer: DGETRS call in LP_BVP_SOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
!mick eff 3/22/2017
        DO Q = 1, N_WEIGHTFUNCS
          NCON(1:NSTREAMS,N,Q) = SCOL2_WF(1:NSTREAMS,Q)
          PCON(1:NSTREAMS,N,Q) = SCOL2_WF(NS1:NS2,Q)
        ENDDO

      ENDIF

!  linearized BVP results
!  ======================

!  Associated quantities

!mick eff 3/22/2017
      DO Q = 1, N_WEIGHTFUNCS
        DO N = 1, NLAYERS
          DO AA = 1, NSTREAMS
            NCON_XVEC(1:NSTREAMS_2,AA,N,Q) = XPOS(1:NSTREAMS_2,AA,N) * NCON(AA,N,Q)
            PCON_XVEC(1:NSTREAMS_2,AA,N,Q) = XNEG(1:NSTREAMS_2,AA,N) * PCON(AA,N,Q)
          ENDDO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE LP_BVP_SOLUTION_MASTER

!

SUBROUTINE LP_BVP_COLUMN_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input, Flags
           DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,    & ! Input, Flags
           DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, TAYLOR_ORDER,         & ! Input, Flag and order
           FOURIER, IBEAM, NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL,        & ! Input, Numbers
           MODIFIED_BCL3, MODIFIED_BCL4, LAYER_TO_VARY, N_LAYER_WFS,     & ! Input
           DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF, QUAD_STRMWTS,   & ! Input, optical and control
           SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SLTERM,       & ! Input, Surface Stuff
           DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! Input, Beam Quantities
           T_DELT_EIGEN, XPOS, XNEG, LCON, MCON, ATERM_SAVE, BTERM_SAVE, & ! Input, Homogeneous/Greens
           CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,           & ! Input, Greens Function
           LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,         & ! Input, Linearized Beam Quantities
           LP_TRANS_ATMOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG,     & ! Input, Linearized Homogeneous
           L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,           & ! Input, Linearized Greens + thermal
           L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                         ! Output

!  Linearized column vector setup (profile weighting functions)
! 4/9/19. Add water-leaving control, reflected direct beam, surface leaving linearization contribution

!    -- 2/28/21. Version 3.8.3. BRDF_F defined locally (remove MAXMOMENTS dimension)
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC arguments to input list

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Surface BRDF and inclusion flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_BRDF_SURFACE
      LOGICAL  , intent(in)  ::  DO_WATER_LEAVING

!  Emission and solar source flags

      LOGICAL  , intent(in)  ::  DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  ::  DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  ::  DO_PLANE_PARALLEL

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  ::  TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS, NTOTAL

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Boundary condition flags

      LOGICAL  , intent(in)  :: MODIFIED_BCL3, MODIFIED_BCL4

!  Linearization control

      INTEGER  , intent(in)  :: LAYER_TO_VARY, N_LAYER_WFS

!  Fourier component, beam number

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Surface stuff
!  -------------

!  Factor, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams
!    2/28/21. Version 3.8.3. BRDF Fourier inputs are defined locally (remove MAXMOMENTS)

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  4/9/19 This is the reflected beam (excludes SLTERM)

      REAL(fpk), intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SLTERM         ( MAXSTREAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Homogeneous
!  -----------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

      REAL(fpk), intent(in)  :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)
!  4/9/19 Add linearization of TRANS_ATMOS for the waterleaving contribution

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_TRANS_ATMOS    ( MAXBEAMS,  MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized Beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COL2_WF    (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  help variables

      INTEGER    :: Q, N, N1, I, I1, CM, C0, K, M
      REAL(fpk)  :: L_HOM, L_PARTIC, LPTERM

!  Local linearized reflectance arrays

      REAL(fpk)  :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
!mick eff 3/22/2017
      INTEGER    :: CNS1, CNS2
      REAL(fpk)  :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)

!  Boundary options

      LOGICAL    :: REGULAR_BCL3, REGULAR_BCL4

!  Initialise
!  ----------

!  Proxies (Layer to vary, Fourier)

      M = FOURIER
      K = LAYER_TO_VARY

!  Zero the results vectors

      COL2_WF(1:NTOTAL,1:MAX_ATMOSWFS) = ZERO

!  Copy already existing thermal linearizations
!    This is a very important zeroing.................!!!!!

      L_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_LAYER_WFS) = ZERO
      L_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_LAYER_WFS) = ZERO
      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO Q = 1, N_LAYER_WFS
          L_WUPPER(1:NSTREAMS_2,K,Q) = L_T_WUPPER(1:NSTREAMS_2,K,Q)
          L_WLOWER(1:NSTREAMS_2,K,Q) = L_T_WLOWER(1:NSTREAMS_2,K,Q)
        ENDDO
      ENDIF

!  Get the linearized beam solution for the first layer
!    -- 2/28/21. Version 3.8.3. Insert CFUNC and DFUNC arguments to input list

      IF ( DO_SOLAR_SOURCES ) THEN
        CALL LP_BEAMSOLUTION_NEQK  &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, K, N_LAYER_WFS, M, IBEAM,       & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
      ENDIF

!  Complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

      IF ( MODIFIED_BCL3 ) THEN

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS
            L_PARTIC = - L_WUPPER(I,N,Q)
            HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
            L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N), L_XPOS(I,1:NSTREAMS,N,Q) ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP1 )
            COL2_WF(I,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

!  No variation case (BCL1)

      ELSE

        COL2_WF(1:NSTREAMS,1:N_LAYER_WFS) = ZERO

      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------------

!  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN

!  .. nothing varying in these layers

        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0  = N1*NSTREAMS_2 - NSTREAMS
          CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS_2
          COL2_WF(CNS1:CNS2,1:N_LAYER_WFS) = ZERO
        ENDDO

      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*NSTREAMS_2 - NSTREAMS

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS_2
            CM = C0 + I
            L_PARTIC = + L_WUPPER(I,N,Q)
            HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
            L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N), L_XPOS(I,1:NSTREAMS,N,Q) ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP1 )
            COL2_WF(CM,Q) = L_PARTIC + L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN

        N  = LAYER_TO_VARY
        N1 = N + 1
        C0 = N*NSTREAMS_2 - NSTREAMS

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, K, N_LAYER_WFS, M, IBEAM,   & ! Input, Numbers
             DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF

!  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS_2
            CM = C0 + I
            L_PARTIC = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I,1:NSTREAMS,N)
            L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), L_XNEG(I,1:NSTREAMS,N,Q) )
            COL2_WF(CM,Q) = L_PARTIC - L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final layer
!  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

          N1 = N + 1
          C0 = N*NSTREAMS_2 - NSTREAMS
!mick eff 3/22/2017
          CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS_2

!  Get the linearized beam solution for the next layer

          IF ( DO_SOLAR_SOURCES ) THEN
            CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N1, K, N_LAYER_WFS, M, IBEAM,   & ! Input, Numbers
             DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
          ENDIF

!  .. contributions from beam solution (direct assign). No homog. variation

          DO Q = 1, N_LAYER_WFS
            COL2_WF(CNS1:CNS2,Q) = L_WUPPER(1:NSTREAMS_2,N1,Q) - L_WLOWER(1:NSTREAMS_2,N,Q)
          ENDDO

!  End layer loop

        ENDDO

!  End BCL5 boundary conditions

      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------------

      N = NLAYERS

!  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN

!  Get the linearized downward-reflected term
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally for each Fourier.

        CALL LP_BVP_SURFACE_SETUP                          &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, M,                     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Compute the solution

        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS
            CM = C0 + I
            I1 = I + NSTREAMS
            L_PARTIC = L_WLOWER(I1,N,Q) - R2_L_PARTIC(I,Q)
            HELP1 =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                  + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                  - R2_L_HOMP(I,1:NSTREAMS,Q)
            HELP2 = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
            L_HOM = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                  + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
            COL2_WF(CM,Q) = - L_PARTIC - L_HOM
          ENDDO
        ENDDO

!  Ordinary BCL6 Component loop
 
      ELSE

!  Get the linearized downward-reflected term
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally for each Fourier.

        CALL LP_BVP_SURFACE_SETUP                          &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, M,                     & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Compute the solution

!mick eff 3/22/2017
        C0 = (N-1)*NSTREAMS_2 + NSTREAMS
        DO Q = 1, N_LAYER_WFS
          DO I = 1, NSTREAMS
            CM = C0 + I
            I1 = I + NSTREAMS              
            COL2_WF(CM,Q) = - L_WLOWER(I1,N,Q) + R2_L_PARTIC(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------
      
! 4/9/19. There was a bug here. Formula was only correct for the reflected beam

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        DO Q = 1, N_LAYER_WFS
          CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS
          COL2_WF(CNS1:CNS2,Q) = COL2_WF(CNS1:CNS2,Q) - RF_DIRECT_BEAM(1:NSTREAMS,IBEAM) &
                      * ( L_DELTAU_VERT(Q,LAYER_TO_VARY)*DELTAU_SLANT(N,LAYER_TO_VARY,IBEAM) )
        ENDDO
      ENDIF

!  4/9/19 Add the linearization due to Adjusted waterleaving term

      IF ( DO_WATER_LEAVING .and. FOURIER .eq. 0 ) THEN
         DO Q = 1, N_LAYER_WFS
            LPTERM = LP_TRANS_ATMOS(IBEAM,LAYER_TO_VARY,Q)
            DO I = 1, NSTREAMS
               CM = C0 + I
               COL2_WF(CM,Q) = COL2_WF(CM,Q) + LPTERM * SLTERM(I)
            ENDDO
         ENDDO
      ENDIF
        
!  debug

!      if ( layer_to_vary.eq.23.and.m.eq.3.and.ipartic.eq.1 ) then
!        write(67,*)layer_to_vary
!        do n = 1, ntotal
!           write(67,*)n,COL2_WF(n,1)
!        enddo
!      endif

!  Copy for the one-layer case

!mick eff 3/22/2017
      IF ( NLAYERS .EQ. 1 ) THEN
        SCOL2_WF(1:NTOTAL,1:N_LAYER_WFS) = COL2_WF(1:NTOTAL,1:N_LAYER_WFS)
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_BVP_COLUMN_SETUP

!

SUBROUTINE LP_BVP_SURFACE_SETUP                            &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, NSTREAMS, & ! Input
            NLAYERS, MODIFIED_BCL4, FOURIER,               & ! Input
            SURFACE_FACTOR, ALBEDO, BRDF_F, N_LAYER_WFS,   & ! Input
            QUAD_STRMWTS, T_DELT_EIGEN, L_T_DELT_EIGEN,    & ! Input
            XPOS, L_XPOS, L_XNEG, L_WLOWER,                & ! Input
            R2_L_PARTIC, R2_L_HOMP, R2_L_HOMM )              ! Output

!  Linearized surface reflectance terms

!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally for each Fourier.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                                MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  input arguments
!  ---------------

!  BRDF flag

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Number of weighting functions

      INTEGER  , intent(in)  :: N_LAYER_WFS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Flag for type of boundary condition

      LOGICAL  , intent(in)  :: MODIFIED_BCL4

!  overall surface flag and surface factor, albedo

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, remove 0:MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_PARTIC(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: FACTOR
      INTEGER    :: AA, Q, N, I, M

!  Initial section
!  ---------------

!  Always zero the result to start

!mick eff 3/22/2017
      R2_L_PARTIC(1:NSTREAMS,1:N_LAYER_WFS) = ZERO
      R2_L_HOMP(1:NSTREAMS,1:NSTREAMS,1:N_LAYER_WFS) = ZERO
      R2_L_HOMM(1:NSTREAMS,1:NSTREAMS,1:N_LAYER_WFS) = ZERO

!  FOurier component

      M = FOURIER

!  Return if no albedo

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Set up Auxiliary arrays
!  -----------------------

      N = NLAYERS

!  Particular integral parts

      DO Q = 1, N_LAYER_WFS
        PV_W(1:NSTREAMS,Q) = L_WLOWER(1:NSTREAMS,N,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MODIFIED_BCL4 ) THEN
        DO Q = 1, N_LAYER_WFS
          DO AA = 1, NSTREAMS
            HV_P(1:NSTREAMS,AA,Q) = QUAD_STRMWTS(1:NSTREAMS) &
                                  * ( L_XPOS(1:NSTREAMS,AA,N,Q) *   T_DELT_EIGEN(AA,N) &
                                      + XPOS(1:NSTREAMS,AA,N)   * L_T_DELT_EIGEN(AA,N,Q) )
            HV_M(1:NSTREAMS,AA,Q) = QUAD_STRMWTS(1:NSTREAMS)*L_XNEG(1:NSTREAMS,AA,N,Q)
          ENDDO
        ENDDO
      ENDIF

!  Integrated Downward reflection (Calculation, Lambertian case)
!  -------------------------------------------------------------

     IF ( .not. DO_BRDF_SURFACE ) THEN

!  Amplitude

        FACTOR = SURFACE_FACTOR * ALBEDO

!  Only a solution for the isotropic part.

        IF ( FOURIER .EQ. 0 ) THEN
          DO Q = 1, N_LAYER_WFS

!  Particular solution

            R2_L_PARTIC(1:NSTREAMS,Q) = FACTOR * SUM( PV_W(1:NSTREAMS,Q) )

!  Homogeneous solutions (only for modified BC)

            IF ( MODIFIED_BCL4 ) THEN
              DO AA = 1, NSTREAMS
                R2_L_HOMP(1:NSTREAMS,AA,Q) = FACTOR * SUM ( HV_P(1:NSTREAMS,AA,Q) )
                R2_L_HOMM(1:NSTREAMS,AA,Q) = FACTOR * SUM ( HV_M(1:NSTREAMS,AA,Q) )
              ENDDO
            ENDIF

!  end parameter loop

          ENDDO
        ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!  ----------------------------------------------------------------

      ELSE

!mick eff 3/22/2017 - moved I loop
        DO Q = 1, N_LAYER_WFS

!  particular solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop M=FOURIER index

          DO I = 1, NSTREAMS
            R2_L_PARTIC(I,Q) = SURFACE_FACTOR * DOT_PRODUCT ( PV_W(1:NSTREAMS,Q), BRDF_F(I,1:NSTREAMS) )
          ENDDO

!  homogeneous solutions
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop M=FOURIER index

          IF ( MODIFIED_BCL4 ) THEN
            DO AA = 1, NSTREAMS
              DO I = 1, NSTREAMS
                R2_L_HOMP(I,AA,Q) = SURFACE_FACTOR * DOT_PRODUCT ( HV_P(1:NSTREAMS,AA,Q), BRDF_F(I,1:NSTREAMS) )
                R2_L_HOMM(I,AA,Q) = SURFACE_FACTOR * DOT_PRODUCT ( HV_M(1:NSTREAMS,AA,Q), BRDF_F(I,1:NSTREAMS) )
              ENDDO
            ENDDO
          ENDIF

!  end parameter and stream loops

        ENDDO
        !ENDDO

!  End BRDF clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_BVP_SURFACE_SETUP

!

SUBROUTINE LP_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, N_WEIGHTFUNCS, M, IB,        & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N
!   Profile variation also in this layer (N = K)

!  2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.
!    -- Necessary because now using different linearization definitions of L_ATERM and L_BTERM

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  ::  M, IB, N

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

      REAL(fpk), intent(in)  :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

! mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Local linearized Green's functioN multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, EPS, DELTA, LAM, MULT
      REAL(fpk)  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_MULT(MAX_ATMOSWFS)
      REAL(fpk)  :: AST, BST, L_AST, L_BST, CST, DST, L_CST, L_DST

!  No linearized particular solution beyond the cutoff layer. ALSO -
!  Nothing if layer is inactive (Does not depend on solution saving)

!  2/28/21. Version 3.8.3. Solutions are already initialized, so exit immediately, no need to zero
!           (Solutions may be initialized to the thermal variables)

!      IF ((DO_SOLUTION_SAVING.AND..NOT.DO_LAYER_SCATTERING(M,N)) &
!             .OR. (N .GT.SOLARBEAM_CUTOFF(IB))) THEN
      IF ( .NOT.DO_LAYER_SCATTERING(M,N) .OR. N.GT.SOLARBEAM_CUTOFF(IB) ) THEN
        RETURN
      ENDIF

!  Green's function solution
!  =========================

!  Set up linearizations of GAMMA constants
!  ----------------------------------------

!  Distinguish two cases:
!  ..(a) quasi-spherical for n > 1
!  ..(b) plane-parallel or QS for n=1

!  Linearizations of optical depth integrations
!  Linearized Green function multipliers

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

      DO AA = 1, NSTREAMS

         ZDEL  = T_DELT_EIGEN(AA,N)
         ZWDEL = ZDEL * WDEL

         AST   = CONST * ATERM_SAVE(AA,N) 
         BST   = CONST * BTERM_SAVE(AA,N) 
         CST   = CONST * CFUNC(AA,N) 
         DST   = CONST * DFUNC(AA,N) 

!  Downwelling, Make allowances for Taylor series

        IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
           EPS   = GAMMA_M(AA,N)
           DELTA = DELTAU_VERT(N)
           LAM   = AVERAGE_SECANT(N,IB)
           DO Q = 1, N_WEIGHTFUNCS
              L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_KEG  = L_KEIGEN(AA,N,Q)
              L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
           ENDDO
        ELSE
           MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
           DO Q = 1, N_WEIGHTFUNCS
              L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
              L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
              L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(AA,N)
           ENDDO
        ENDIF

!    -- 2/28/21. Changed code to use regular L_ATERM derivative

        DO Q = 1, N_WEIGHTFUNCS
!           L_AST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_ATERM_SAVE(AA,N,Q)
!           L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
           L_CST =  LP_INITIAL_TRANS(N,N,IB,Q) * CFUNC(AA,N) + L_MULT(Q)
           L_GFUNC_DN(AA,Q) = AST * L_CST + L_ATERM_SAVE(AA,N,Q) * CST
        ENDDO

!  Upwelling

        MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
        DO Q = 1, N_WEIGHTFUNCS
           L_ZDEL =  L_T_DELT_EIGEN  (AA,N,Q)   ; L_KEG  = L_KEIGEN(AA,N,Q)
           L_WDEL =  LP_T_DELT_MUBAR (N,N,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,N,IB,Q)
           L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(AA,N)
        ENDDO

!    -- 2/28/21. Changed code to use regular L_BTERM derivative

        DO Q = 1, N_WEIGHTFUNCS
!           L_BST =  LP_INITIAL_TRANS(N,N,IB,Q)  + L_BTERM_SAVE(AA,N,Q)
!           L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
           L_DST =  LP_INITIAL_TRANS(N,N,IB,Q) * DFUNC(AA,N) + L_MULT(Q)
           L_GFUNC_UP(AA,Q) = BST * L_DST + L_BTERM_SAVE(AA,N,Q) * DST
        ENDDO

!  End discrete ordinate loop

      ENDDO

!  Set linearized form of particular integral at boundaries
!  --------------------------------------------------------

!mick eff 3/22/2017
      DO Q = 1, N_WEIGHTFUNCS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          LBSOL(I,Q,1)  = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q),  XPOS(I1,1:NSTREAMS,N)   ) &
                        + DOT_PRODUCT (   GFUNC_UP(1:NSTREAMS,N),L_XPOS(I1,1:NSTREAMS,N,Q) )
          LBSOL(I1,Q,1) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q),  XPOS(I, 1:NSTREAMS,N)   ) &
                        + DOT_PRODUCT (   GFUNC_UP(1:NSTREAMS,N),L_XPOS(I, 1:NSTREAMS,N,Q) )
          LBSOL(I1,Q,2) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q),  XPOS(I1,1:NSTREAMS,N)   ) &
                        + DOT_PRODUCT (   GFUNC_DN(1:NSTREAMS,N),L_XPOS(I1,1:NSTREAMS,N,Q) )
          LBSOL(I,Q,2)  = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q),  XPOS(I, 1:NSTREAMS,N)   ) &
                        + DOT_PRODUCT (   GFUNC_DN(1:NSTREAMS,N),L_XPOS(I, 1:NSTREAMS,N,Q) )
        ENDDO
      ENDDO

!mick fix 2/17/12 - code added to include both solar & thermal sources.
!                   like LC_BEAMSOLUTION_NEQK now
! Add to existing solution

!mick eff 3/22/2017
      DO Q = 1, N_WEIGHTFUNCS
        L_WUPPER(1:NSTREAMS_2,N,Q) = L_WUPPER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,1)
        L_WLOWER(1:NSTREAMS_2,N,Q) = L_WLOWER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,2)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BEAMSOLUTION_NEQK

!

SUBROUTINE LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER, & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, K, K_PARAMETERS, M, IB,      & ! Input, Numbers
             DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N
!   Profile variation Not in this layer (N =/ K)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2

!  Fourier number, beam index, layer index

      INTEGER  , intent(in)  :: M, IB, N

!  Varying lyaer, number of varying parameters (input)

      INTEGER  , intent(in)  :: K, K_PARAMETERS

!  Input optical properties after delta-M scaling. 

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Local linearized Green's function multipliers

      REAL(fpk)  :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: AA, I, I1, Q
      REAL(fpk)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT, T1
      REAL(fpk)  :: L_WDEL, L_LAM, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)
      REAL(fpk)  :: LBSOL(MAXSTREAMS_2,MAX_ATMOSWFS,2)

!  Green's function solution
!  =========================

!  No linearized particular solution beyond the cutoff layer. ALSO--
!  Nothing if layer is inactive (Does not depend on solution saving)
!  2/28/21. Solutions are already initialized, so exit immediately, no need to initialize

      IF ( .NOT.DO_LAYER_SCATTERING(M,N) .OR. N.GT.SOLARBEAM_CUTOFF(IB) )  RETURN

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

!  Set up linearizations of GAMMA constants

!  Mick note 9/25/2013 - LP_GAMMA_P & LP_GAMMA_M were previously
!  normalized by GAMMA_P & GAMMA_M, respectively.  With the new
!  definitions of GAMMA_P & GAMMA_M (now straight sums and differences
!  instead of their corresponding reciprocals), we define LP_GAMMA_M as
!  follows:

        CONST   = INITIAL_TRANS(N,IB)
        WDEL    = T_DELT_MUBAR(N,IB)

!  Start discrete ordinate loop

        DO AA = 1, NSTREAMS

          ZDEL  = T_DELT_EIGEN(AA,N)
          ZWDEL = ZDEL * WDEL
          AST   = CONST * ATERM_SAVE(AA,N) 
          BST   = CONST * BTERM_SAVE(AA,N) 

!mick fix 9/4/2013 - small numbers analysis added

!  Downwelling, Make allowances for Taylor series

          IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
            EPS   = GAMMA_M(AA,N)
            DELTA = DELTAU_VERT(N)
            LAM   = AVERAGE_SECANT(N,IB)
            DO Q = 1, K_PARAMETERS
              L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, ZERO, ZERO, L_LAM, WDEL, LAM, L_MULT(Q) )
            ENDDO
          ELSE
            MULT = ( ZDEL - WDEL ) / GAMMA_M(AA,N)
            DO Q = 1, K_PARAMETERS
              L_WDEL =  LP_T_DELT_MUBAR (N,K,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
              L_MULT(Q) = ( - L_WDEL - MULT * L_LAM ) / GAMMA_M(AA,N)
            ENDDO
          ENDIF
          DO Q = 1, K_PARAMETERS
            L_AST =  LP_INITIAL_TRANS(N,K,IB,Q) 
            L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * L_AST + L_MULT(Q) * AST
          ENDDO

!  Upwelling

          MULT = ( ONE - ZWDEL ) / GAMMA_P(AA,N)
          DO Q = 1, K_PARAMETERS
            L_WDEL =  LP_T_DELT_MUBAR (N,K,IB,Q) ; L_LAM  = LP_AVERAGE_SECANT(N,K,IB,Q)
            L_MULT(Q) = - ( L_WDEL*ZDEL + MULT*L_LAM ) / GAMMA_P(AA,N)
          ENDDO
          DO Q = 1, K_PARAMETERS
            L_BST =  LP_INITIAL_TRANS(N,K,IB,Q)
            L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * L_BST + L_MULT(Q) * BST
          ENDDO

!  End eigenloop

        ENDDO

!  ..(b) plane-parallel and qs for n = 1

      ELSE

!mick eff 3/22/2017
        DO Q = 1, K_PARAMETERS
          T1 = LP_INITIAL_TRANS(N,K,IB,Q)
          L_GFUNC_DN(1:NSTREAMS,Q) = GFUNC_DN(1:NSTREAMS,N) * T1
          L_GFUNC_UP(1:NSTREAMS,Q) = GFUNC_UP(1:NSTREAMS,N) * T1
        ENDDO

      ENDIF

!  Set linearized form of particular integral

!mick eff 3/22/2017
      DO Q = 1, K_PARAMETERS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          LBSOL(I,Q,1)  = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q), XPOS(I1,1:NSTREAMS,N) )
          LBSOL(I1,Q,1) = DOT_PRODUCT ( L_GFUNC_UP(1:NSTREAMS,Q), XPOS(I, 1:NSTREAMS,N) )
          LBSOL(I1,Q,2) = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q), XPOS(I1,1:NSTREAMS,N) )
          LBSOL(I,Q,2)  = DOT_PRODUCT ( L_GFUNC_DN(1:NSTREAMS,Q), XPOS(I, 1:NSTREAMS,N) )
        ENDDO
      ENDDO

!mick fix 2/17/12 - code added to include both solar & thermal sources.
! Add to existing solution

!mick eff 3/22/2017
      DO Q = 1, K_PARAMETERS
        L_WUPPER(1:NSTREAMS_2,N,Q) = L_WUPPER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,1)
        L_WLOWER(1:NSTREAMS_2,N,Q) = L_WLOWER(1:NSTREAMS_2,N,Q) + LBSOL(1:NSTREAMS_2,Q,2)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BEAMSOLUTION_NNEK

!

SUBROUTINE LP_BVPTEL_SOLUTION_MASTER &
       ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,         & ! Input, Flags
         DO_PLANE_PARALLEL, TAYLOR_ORDER, FOURIER, IBEAM, NLAYERS, NSTREAMS,     & ! Input, Numbers
         NSTREAMS_2, LAYER_TO_VARY, N_LAYER_WFS, N_SUPDIAG, N_SUBDIAG,           & ! Input, Numbers
         ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,                         & ! Input, Numbers
         DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                           & ! Input, optical and control
         SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT,        & ! Input, Surface +DB Stuff
         INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_DELT_DISORDS,            & ! Input, Beam Quantities
         T_DELT_EIGEN, XPOS, XNEG, WLOWER, LCON, MCON, ATERM_SAVE, BTERM_SAVE,   & ! Input, Homogeneous + PI
         CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                     & ! Input, Greens Function
         BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, LCON_XVEC, MCON_XVEC,           & ! Input, BVP Bandmat
         LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, L_T_DELT_DISORDS, & ! Input , Linearized Beam Quantities
         L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input , Linearized solutions
         NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,                   & ! Output, Linearized Constants + PI
         STATUS, MESSAGE, TRACE )                                                  ! Output - Exception handling

!  Linearization of the Telescoped Boundary Problem Solution
!   Version 3.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, each Fourier component
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL,  &
                                ZERO, LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  surface control. Always going to be BRDF here

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: layer varying and number of varying parameters (input)

      INTEGER  , intent(in)  ::  LAYER_TO_VARY, N_LAYER_WFS

!  BVProblem Band matrix control

      INTEGER  , intent(in)  ::  N_SUBDIAG, N_SUPDIAG

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Surface Quantities
!  ------------------

!  surface control

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Fourier components of BRDF
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop 0:MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  BVProblem inputs
!  ----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2       (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDTELMAT2 (MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT    (MAXSTREAMS_2)

!  discrete ordinate factors (BVP telescoping, solutions saving)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Particular integral at lower boundary

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

      REAL(fpk), intent(in)  :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(out) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving linearized BCs

      REAL(fpk)  :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk)  :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local help variables 

      INTEGER    :: I, I1, K, Q, M
      INTEGER    :: NS, N, N1, NAF, NAL, NAL1, AA, C0
      REAL(fpk)  :: SHOM
!mick eff 3/22/2017
      INTEGER    :: NS1, NS2, CM1, CM2, CP1, CP2
      REAL(fpk)  :: HELP1 (NSTREAMS)

!  Initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Wf count

      M = FOURIER
      Q = 0
!mick eff 3/22/2017
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Set up linearized BVP, column vector  B' where AX = B'
!  ======================================================

!  Bulk: Compute the main column B' where AX = B'

!   2/28/21. Version 3.8.3. BRDF_F array defined locally, for each Fourier component
!   2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs. Rearrange I/O slightly.

      CALL LP_BVPTEL_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,  & ! Input, Flags and order
             DO_PLANE_PARALLEL, TAYLOR_ORDER, NLAYERS, NSTREAMS, NSTREAMS_2,  & ! Input, Flags and order
             LAYER_TO_VARY, N_LAYER_WFS, FOURIER, IBEAM,                      & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,                  & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                    & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Input, Beam Quantities
             T_DELT_EIGEN, T_DELT_DISORDS, XPOS, XNEG, WLOWER, LCON, MCON,    & ! Input, Homogeneous + WLOWER
             SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT, & ! Input, Surface + Direct-beam
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,              & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,              & ! Input, Linearized Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,            & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, L_T_DELT_DISORDS,      & ! Input, Linearized Homogeneous 
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                         ! Output

!  Solve linearized BVP, several active layers
!  ===========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  BV solution for linearized integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, N_LAYER_WFS, &
               BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2_WF, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = 'DGBTRS call in L_BVPTEL_BEAMSOLUTION_MASTER'
         STATUS  = LIDORT_SERIOUS
         RETURN
        ENDIF

!  Set linearized integration constants, active layers

!mick eff 3/22/2017
        C0 = - NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          CM1 = C0 + 1         ; CM2 = C0 + NSTREAMS
          CP1 = CM1 + NSTREAMS ; CP2 = CM2 + NSTREAMS
          DO Q = 1, N_LAYER_WFS
            NCON(1:NSTREAMS,N,Q) = COLTEL2_WF(CM1:CM2,Q)
            PCON(1:NSTREAMS,N,Q) = COLTEL2_WF(CP1:CP2,Q)
          ENDDO
        ENDDO

!  Solve linearized BVP: Single Layer only
!  =======================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ('N', NSTREAMS_2, N_LAYER_WFS,  & 
                    SMAT2, MAXSTREAMS_2, SIPIVOT, SCOL2_WF, MAXSTREAMS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in LP_BVPTEL_BEAMSOLUTION_MASTER'
          STATUS  = LIDORT_SERIOUS
        ENDIF

!  Set linearized integration constants for active layer

!mick eff 3/22/2017
        N = ACTIVE_LAYERS(1)
        DO Q = 1, N_LAYER_WFS
          NCON(1:NSTREAMS,N,Q) = SCOL2_WF(1:NSTREAMS,Q)
          PCON(1:NSTREAMS,N,Q) = SCOL2_WF(NS1:NS2,Q)
        ENDDO

!  End clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

!mick eff 3/22/2017
      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO Q = 1, N_LAYER_WFS
          DO K = 1, NSTREAMS
            NCON_XVEC(1:NSTREAMS_2,K,N,Q) = NCON(K,N,Q) * XPOS(1:NSTREAMS_2,K,N)
            PCON_XVEC(1:NSTREAMS_2,K,N,Q) = PCON(K,N,Q) * XNEG(1:NSTREAMS_2,K,N)
          ENDDO
        ENDDO
      ENDDO

!  Set linearized integration constants for non-active layers
!  ==========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required if the first active
!       layer is the varying layer


      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        IF ( LAYER_TO_VARY.EQ.NAF ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              HELP1 = PCON_XVEC(I1,1:NSTREAMS,NAF,Q) + MCON(1:NSTREAMS,NAF)*L_XNEG(I1,1:NSTREAMS,NAF,Q)
              SHOM  = SUM ( NCON_XVEC(I1,1:NSTREAMS,NAF,Q) ) &
                    + DOT_PRODUCT ( LCON(1:NSTREAMS,NAF), L_XPOS(I1,1:NSTREAMS,NAF,Q) ) &
                    + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAF), HELP1 ) &
                    + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAF,Q), MCON_XVEC(I1,1:NSTREAMS,NAF) )
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY.LT.NAF) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              SHOM  = SUM ( NCON_XVEC(I1,1:NSTREAMS,NAF,Q) ) &
                    + DOT_PRODUCT ( T_DELT_EIGEN(1:NSTREAMS,NAF), PCON_XVEC(I1,1:NSTREAMS,NAF,Q) )
              PCON(I,N1,Q) = L_WUPPER(I1,NAF,Q) + SHOM
              NCON(I,N1,Q) = ZERO
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            PCON(1:NSTREAMS,N1,Q) = ZERO
            NCON(1:NSTREAMS,N1,Q) = ZERO
          ENDDO
        ENDIF
      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards.
!   Additional linearizations if you are passing through the varying layer.

!mick eff 3/22/2017
      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO Q = 1, N_LAYER_WFS
          NCON(1:NSTREAMS,N,Q) = ZERO
          PCON(1:NSTREAMS,N,Q) = T_DELT_DISORDS(1:NSTREAMS,N1) * PCON(1:NSTREAMS,N1,Q)
          IF ( N1 .EQ. LAYER_TO_VARY ) THEN
            PCON(1:NSTREAMS,N,Q) = PCON(1:NSTREAMS,N,Q) + L_T_DELT_DISORDS(1:NSTREAMS,N1,Q) * MCON(1:NSTREAMS,N1)
          ENDIF
        ENDDO
      ENDDO     

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-constants, always required

        IF ( LAYER_TO_VARY .EQ. NAL ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              HELP1 = NCON_XVEC(I,1:NSTREAMS,NAL,Q) + LCON(1:NSTREAMS,NAL)*L_XPOS(I,1:NSTREAMS,NAL,Q)
              SHOM  = SUM ( PCON_XVEC(I,1:NSTREAMS,NAL,Q) ) &
                    + DOT_PRODUCT ( MCON(1:NSTREAMS,NAL), L_XNEG(I,1:NSTREAMS,NAL,Q) ) &
                    + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 ) &
                    + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAL,Q), LCON_XVEC(I,1:NSTREAMS,NAL) )
              NCON(I,NAL1,Q) = L_WLOWER(I,NAL,Q) + SHOM
            ENDDO
          ENDDO
        ELSE IF ( LAYER_TO_VARY .LT. NAL ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              HELP1 = NCON_XVEC(I,1:NSTREAMS,NAL,Q)
              SHOM  = SUM ( PCON_XVEC(I,1:NSTREAMS,NAL,Q) ) &
                    + DOT_PRODUCT ( T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 )
              NCON(I,NAL1,Q) = L_WLOWER(I,NAL,Q) + SHOM
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              HELP1 = NCON_XVEC(I,1:NSTREAMS,NAL,Q)
              SHOM  = SUM ( PCON_XVEC(I,1:NSTREAMS,NAL,Q) ) &
                    + DOT_PRODUCT ( T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 )
              NCON(I,NAL1,Q) = SHOM
            ENDDO
          ENDDO
        ENDIF

!mick eff 3/22/2017
        DO N = NAL + 2, NLAYERS
          N1 = N - 1
          DO Q = 1, N_LAYER_WFS
            NCON(1:NSTREAMS,N,Q) = T_DELT_DISORDS(1:NSTREAMS,N1) * NCON(1:NSTREAMS,N1,Q)
            IF ( N1 .EQ. LAYER_TO_VARY ) THEN
              NCON(1:NSTREAMS,N,Q) = NCON(1:NSTREAMS,N,Q) + L_T_DELT_DISORDS(1:NSTREAMS,N1,Q) *  LCON(1:NSTREAMS,N1)
            ENDIF
          ENDDO
        ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.

        IF ( DO_INCLUDE_SURFACE ) THEN
          IF ( LAYER_TO_VARY .EQ. NAL ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                HELP1 = NCON_XVEC(I1,1:NSTREAMS,NAL,Q) +  LCON(1:NSTREAMS,NAL)*L_XPOS(I1,1:NSTREAMS,NAL,Q)
                SHOM  = SUM ( PCON_XVEC(I1,1:NSTREAMS,NAL,Q) ) &
                      + DOT_PRODUCT ( MCON(1:NSTREAMS,NAL), L_XNEG(I1,1:NSTREAMS,NAL,Q) ) &
                      + DOT_PRODUCT (   T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 ) &
                      + DOT_PRODUCT ( L_T_DELT_EIGEN(1:NSTREAMS,NAL,Q), LCON_XVEC(I1,1:NSTREAMS,NAL)  )
                PCON(I,NAL1,Q) = ( L_WLOWER(I1,NAL,Q) + SHOM ) / T_DELT_DISORDS(I,NAL1)
              ENDDO
            ENDDO
          ELSE IF ( LAYER_TO_VARY .LT. NAL ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                HELP1 = NCON_XVEC(I1,1:NSTREAMS,NAL,Q)
                SHOM  = SUM ( PCON_XVEC(I1,1:NSTREAMS,NAL,Q) ) &
                      + DOT_PRODUCT ( T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 )
                PCON(I,NAL1,Q) = ( L_WLOWER(I1,NAL,Q) + SHOM ) / T_DELT_DISORDS(I,NAL1)
              ENDDO
            ENDDO
          ELSE
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                HELP1 = NCON_XVEC(I1,1:NSTREAMS,NAL,Q)
                SHOM  = SUM ( PCON_XVEC(I1,1:NSTREAMS,NAL,Q) ) &
                      + DOT_PRODUCT ( T_DELT_EIGEN(1:NSTREAMS,NAL), HELP1 )
                IF ( LAYER_TO_VARY.EQ.NAL1 ) THEN
                  PCON(I,NAL1,Q) = ( SHOM - L_T_DELT_DISORDS(I,NAL1,Q) * MCON(I,NAL1) ) / T_DELT_DISORDS(I,NAL1)
                ELSE
                  PCON(I,NAL1,Q) = SHOM / T_DELT_DISORDS(I,NAL1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          DO N = NAL + 2, NLAYERS
            N1 = N - 1
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              IF ( N.EQ.LAYER_TO_VARY ) THEN
                PCON(1:NSTREAMS,N,Q) = ( PCON(1:NSTREAMS,N1,Q) &
                                     - L_T_DELT_DISORDS(1:NSTREAMS,N,Q) * MCON(1:NSTREAMS,N) ) / T_DELT_DISORDS(1:NSTREAMS,N)
              ELSE
                PCON(1:NSTREAMS,N,Q) = PCON(1:NSTREAMS,N1,Q) / T_DELT_DISORDS(1:NSTREAMS,N)
              ENDIF
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          PCON(1:NSTREAMS,NAL1:NLAYERS,1:N_LAYER_WFS) = ZERO
        ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  Associated quantities for inactive layers
!  -----------------------------------------

!  Atmosphere layers with no scattering

!mick eff 3/22/2017
      DO Q = 1, N_LAYER_WFS
        DO N = 1, NLAYERS
          IF ( N.LT.NAF .OR. N.GT.NAL ) THEN
            DO AA = 1, NSTREAMS
              NCON_XVEC(1:NSTREAMS_2,AA,N,Q) = NCON(AA,N,Q) * XPOS(1:NSTREAMS_2,AA,N)
              PCON_XVEC(1:NSTREAMS_2,AA,N,Q) = PCON(AA,N,Q) * XNEG(1:NSTREAMS_2,AA,N)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BVPTEL_SOLUTION_MASTER

!

SUBROUTINE LP_BVPTEL_SURFACE_SETUP                                             &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, FOURIER,                 & ! Input
            LAYER_TO_VARY, N_LAYER_WFS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,  & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  Linearized surface reflectance terms
!    Version 3.8. Major extension to BRDFs with TELESCOPING. Implemented May 2016.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                                MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Surface (BRDF) flag

      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS, N

!  Layer to vary and Number of weighting functions

      INTEGER  , intent(in)  :: LAYER_TO_VARY
      INTEGER  , intent(in)  :: N_LAYER_WFS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop 0:MAXMOMENTS dimensioning 

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  transmittance factors discrete ordinate streams, and Linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output arguments
!  ----------------

      REAL(fpk), intent(out) :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface
!     Calculated here, but  could be done earlier and passed in

      REAL(fpk), intent(out) :: CUMTRANS(MAXSTREAMS)
      REAL(fpk), intent(out) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      REAL(fpk)  :: QCUMTRANS(MAXSTREAMS)
      REAL(fpk)  :: L_QCUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: PS_W(MAXSTREAMS), PV_W(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HS_P(MAXSTREAMS,MAXSTREAMS), HV_P(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: HS_M(MAXSTREAMS,MAXSTREAMS), HV_M(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: BEAM, L_BEAM, HOMP, HOMM, L_HOMP, L_HOMM
      INTEGER    :: AA, Q, I, M, N1
      LOGICAL    :: MBOUNDARY
!mick eff 3/22/2017
      REAL(fpk)  :: HSP(NSTREAMS), HSM(NSTREAMS), L_HSP(NSTREAMS), L_HSM(NSTREAMS)

!  Initial section
!  ---------------

!  Always zero the result to start

!mick eff 3/22/2017

      R2_L_BEAM(1:NSTREAMS,1:N_LAYER_WFS) = ZERO
      R2_L_HOMP(1:NSTREAMS,1:NSTREAMS,1:N_LAYER_WFS) = ZERO
      R2_L_HOMM(1:NSTREAMS,1:NSTREAMS,1:N_LAYER_WFS) = ZERO

      CUMTRANS(1:NSTREAMS) = ONE
      L_CUMTRANS(1:NSTREAMS,1:N_LAYER_WFS) = ZERO

!  Proxy Fourier component

      M = FOURIER

!  Boundary flag

      MBOUNDARY = .true.

!  Return if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Cumulative transmittance

!mick eff 3/22/2017 - initializing turned off (already done above)
!                   - IF block
      !CUMTRANS(1:NSTREAMS) = ONE ; L_CUMTRANS(1:NSTREAMS,:) = ZERO
      DO N1 = NLAYERS, N+1, -1
         DO Q = 1, N_LAYER_WFS
           IF ( N1.eq.LAYER_TO_VARY ) THEN
             L_CUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * T_DELT_DISORDS(1:NSTREAMS,N1) &
                                      +   CUMTRANS(1:NSTREAMS) * L_T_DELT_DISORDS(1:NSTREAMS,N1,Q)
           ELSE
             L_CUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * T_DELT_DISORDS(1:NSTREAMS,N1)
           ENDIF
         ENDDO
         CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
      ENDDO

!  Stored

      QCUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * QUAD_STRMWTS(1:NSTREAMS)
      DO Q = 1, N_LAYER_WFS
         L_QCUMTRANS(1:NSTREAMS,Q) = L_CUMTRANS(1:NSTREAMS,Q) * QUAD_STRMWTS(1:NSTREAMS)
      ENDDO

!  Particular integral

!mick eff 3/22/2017
      PS_W(1:NSTREAMS) = WLOWER(1:NSTREAMS,N) * QCUMTRANS(1:NSTREAMS)
      DO Q = 1, N_LAYER_WFS
        PV_W(1:NSTREAMS,Q) = L_WLOWER(1:NSTREAMS,N,Q) *   QCUMTRANS(1:NSTREAMS) &
                           +   WLOWER(1:NSTREAMS,N)   * L_QCUMTRANS(1:NSTREAMS,Q)
      ENDDO

!    Modified boundary condition: homogeneous parts

      IF ( MBOUNDARY ) THEN
        IF ( N.EQ.LAYER_TO_VARY ) THEN
!mick eff 3/22/2017
          DO AA = 1, NSTREAMS
            HSP = XPOS(1:NSTREAMS,AA,N) * T_DELT_EIGEN(AA,N)
            HSM = XNEG(1:NSTREAMS,AA,N)
            HS_P(1:NSTREAMS,AA) = HSP * QCUMTRANS(1:NSTREAMS)
            HS_M(1:NSTREAMS,AA) = HSM * QCUMTRANS(1:NSTREAMS)
            DO Q = 1, N_LAYER_WFS
              L_HSP = L_XPOS(1:NSTREAMS,AA,N,Q) * T_DELT_EIGEN(AA,N) + XPOS(1:NSTREAMS,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
              L_HSM = L_XNEG(1:NSTREAMS,AA,N,Q)
              HV_P(1:NSTREAMS,AA,Q) = L_HSP * QCUMTRANS(1:NSTREAMS) + HSP * L_QCUMTRANS(1:NSTREAMS,Q)
              HV_M(1:NSTREAMS,AA,Q) = L_HSM * QCUMTRANS(1:NSTREAMS) + HSM * L_QCUMTRANS(1:NSTREAMS,Q)
            ENDDO
          ENDDO
        ELSE
!mick eff 3/22/2017
          DO AA = 1, NSTREAMS
            HSP = XPOS(1:NSTREAMS,AA,N) * T_DELT_EIGEN(AA,N)
            HSM = XNEG(1:NSTREAMS,AA,N)
            HS_P(1:NSTREAMS,AA) = HSP * QCUMTRANS(1:NSTREAMS)
            HS_M(1:NSTREAMS,AA) = HSM * QCUMTRANS(1:NSTREAMS)
            DO Q = 1, N_LAYER_WFS
              HV_P(1:NSTREAMS,AA,Q) = HSP * L_QCUMTRANS(1:NSTREAMS,Q)
              HV_M(1:NSTREAMS,AA,Q) = HSM * L_QCUMTRANS(1:NSTREAMS,Q)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Integrated Downward reflection (Calculation, Bidirectional case)
!     homogeneous and particular solutions.
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop M=FOURIER index 

      DO I = 1, NSTREAMS
        BEAM = DOT_PRODUCT( PS_W(1:NSTREAMS),BRDF_F(I,1:NSTREAMS) )
        DO Q = 1, N_LAYER_WFS
          L_BEAM = DOT_PRODUCT( PV_W(1:NSTREAMS,Q),BRDF_F(I,1:NSTREAMS) )
          R2_L_BEAM(I,Q) = SURFACE_FACTOR * ( L_BEAM * CUMTRANS(I) + BEAM * L_CUMTRANS(I,Q) )
        ENDDO
        IF ( MBOUNDARY ) THEN
          DO AA = 1, NSTREAMS
            HOMP = DOT_PRODUCT(HS_P(1:NSTREAMS,AA),BRDF_F(I,1:NSTREAMS) )
            HOMM = DOT_PRODUCT(HS_M(1:NSTREAMS,AA),BRDF_F(I,1:NSTREAMS) )
            DO Q = 1, N_LAYER_WFS
              L_HOMP = DOT_PRODUCT(HV_P(1:NSTREAMS,AA,Q),BRDF_F(I,1:NSTREAMS) )
              L_HOMM = DOT_PRODUCT(HV_M(1:NSTREAMS,AA,Q),BRDF_F(I,1:NSTREAMS) )
              R2_L_HOMP(I,AA,Q) = SURFACE_FACTOR * ( L_HOMP * CUMTRANS(I) + HOMP * L_CUMTRANS(I,Q) )
              R2_L_HOMM(I,AA,Q) = SURFACE_FACTOR * ( L_HOMM * CUMTRANS(I) + HOMM * L_CUMTRANS(I,Q) )
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_BVPTEL_SURFACE_SETUP

!

SUBROUTINE LP_BVPTEL_COLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_LAYER_SCATTERING,  & ! Input, Flags and order
             DO_PLANE_PARALLEL, TAYLOR_ORDER, NLAYERS, NSTREAMS, NSTREAMS_2,  & ! Input, Numbers
             LAYER_TO_VARY, N_LAYER_WFS, FOURIER, IBEAM,                      & ! Input, Numbers
             ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,                  & ! Input, Numbers for Telescoping
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,                    & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Input, Beam Quantities
             T_DELT_EIGEN, T_DELT_DISORDS, XPOS, XNEG, WLOWER, LCON, MCON,    & ! Input, Homogeneous + WLOWER
             SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM, DELTAU_SLANT, & ! Input, Surface + Direct-beam
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,              & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,              & ! Input, Linearized Greens Function
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,            & ! Input, Linearized Beam Quantities
             L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, L_T_DELT_DISORDS,      & ! Input, Linearized Homogeneous
             L_WUPPER, L_WLOWER, COLTEL2_WF, SCOL2_WF )                         ! Output

!  Column setup for the linearized telescoped BVP

!  2/28/21. Version 3.8.3. BRDF_F array defined locally, for each Fourier component 
!  2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs. Rearrange I/O slightly.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  input arguments
!  ===============

!  Control and Optical
!  -------------------

!  Surface (BRDF) and directbeam flag

      LOGICAL  , intent(in)  ::  DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  ::  DO_INCLUDE_DIRECTBEAM

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  ::  DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of layers

      INTEGER  , intent(in)  ::  NLAYERS

!  Number of streams

      INTEGER  , intent(in)  ::  NSTREAMS, NSTREAMS_2

!  Linearization control: layer and number of varying parameters (input)

      INTEGER  , intent(in)  ::  LAYER_TO_VARY, N_LAYER_WFS

!  Fourier number, beam index

      INTEGER  , intent(in)  ::  FOURIER, IBEAM

!  Number of telescoped layers

      INTEGER  , intent(in)  ::  NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  ::  ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  ::  N_BVTELMATRIX_SIZE

!  Input optical properties after delta-M scaling. 
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Surface Quantities
!  ------------------

!  surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, drop 0:MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Quadrature input

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Direct beam solutions

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Solution variables
!  ------------------

!  transmittance factors discrete ordinate streams, and Linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  particular solutions at Lower boundaries

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  Green functions
!  ---------------

!  Green function Multipliers for solution
!    -- 2/28/21. Version 3.8.3. Introduce CFUNC and DFUNC inputs.

      REAL(fpk), intent(in)  :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Linearized
!  ----------

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LP_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvalues, Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(fpk), intent(inout) :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving linearized BCs

      REAL(fpk), intent(out) :: COLTEL2_WF (MAXTOTAL,    MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: SCOL2_WF   (MAXSTREAMS_2,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Output arguments from the Surface setup (reflectances and cumulative transmittances)
!  cumulative tranmsittance (and linearization) from bottom of lowest active layer to surface

      REAL(fpk) :: R2_L_BEAM(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: R2_L_HOMP(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: R2_L_HOMM(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

      REAL(fpk) :: CUMTRANS(MAXSTREAMS)
      REAL(fpk) :: L_CUMTRANS(MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER    :: Q,N,N1,NS,I,I1,CM,C0,NAF,M
      REAL(fpk)  :: BEAM, L_HOM, L_BEAM
!mick eff 3/22/2017
      INTEGER     :: CNS1, CNS2
      REAL(fpk)   :: HELP1 (NSTREAMS), HELP2 (NSTREAMS)

!  Try this safety-first zeroing

!mick eff 3/22/2017
      L_WUPPER(1:NSTREAMS_2,1:NLAYERS,1:N_LAYER_WFS) = ZERO
      L_WLOWER(1:NSTREAMS_2,1:NLAYERS,1:N_LAYER_WFS) = ZERO

!  Fourier proxy

      M = FOURIER

!  Get the linearized solutions for all active layers
!    Always need this, regardless of number of active layers

!   2/28/21. Version 3.8.3. Insert CFUNC and DFUNC arguments to input list

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        IF ( N.EQ.LAYER_TO_VARY ) THEN
          CALL LP_BEAMSOLUTION_NEQK &
           ( DO_LAYER_SCATTERING, TAYLOR_ORDER,                    & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, N_LAYER_WFS, M, IBEAM,       & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_CUTOFF,         & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, & ! Input, Homogeneous solution stuff
             CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output
        ENDIF
      ENDDO

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   GENERAL CASE, 2 or MORE LAYERS
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  Zero column vector

!mick eff 3/22/2017
        !COLTEL2_WF(1:N_BVTELMATRIX_SIZE,:) = ZERO
        COLTEL2_WF(1:N_BVTELMATRIX_SIZE,1:N_LAYER_WFS) = ZERO

!  Top of first active layer, first boundary condition
!  ---------------------------------------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)
        NAF = N

!  If this active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

        IF ( LAYER_TO_VARY .EQ. N ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              L_BEAM = - L_WUPPER(I,N,Q)
              HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N),L_XPOS(I,1:NSTREAMS,N,Q) ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N),HELP1 )
              COLTEL2_WF(I,Q) = L_BEAM - L_HOM
            ENDDO
          ENDDO

!  Otherwise if varying layer is above first active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,  & ! Input, Numbers
             DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            COLTEL2_WF(1:NSTREAMS,Q) = - L_WUPPER(1:NSTREAMS,N,Q)
          ENDDO

        ENDIF

!  Intermediate boundaries between active layers
!  ---------------------------------------------

        DO NS = 1, NLAYERS_TEL - 1

!  Offsets

          N  = ACTIVE_LAYERS(NS)
          N1 = N + 1
          C0 = NS*NSTREAMS_2 - NSTREAMS

!  if N is the varying layer, immediately above boundary
!  Get the linearized beam solution for the next layer N1

          IF ( N .EQ. LAYER_TO_VARY ) THEN

            CALL LP_BEAMSOLUTION_NNEK &
              ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
                NSTREAMS, NSTREAMS_2, N1, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM, & ! Input, Numbers
                DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
                INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
                LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
                T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
                GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
                L_WUPPER, L_WLOWER )                                    ! Output

!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
                HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I,1:NSTREAMS,N,Q) &
                         + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I,1:NSTREAMS,N)
                L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                         + DOT_PRODUCT ( MCON(1:NSTREAMS,N), L_XNEG(I,1:NSTREAMS,N,Q) )
                COLTEL2_WF(CM,Q) = L_BEAM - L_HOM
              ENDDO
            ENDDO

!  If N1 is the varying layer, immediately below boundary
!    Only require contributions from this layer

          ELSE IF ( N1 .EQ. LAYER_TO_VARY ) THEN

!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                L_BEAM = + L_WUPPER(I,N1,Q)
                HELP1    =   T_DELT_EIGEN(1:NSTREAMS,N1)   * L_XNEG(I,1:NSTREAMS,N1,Q) &
                         + L_T_DELT_EIGEN(1:NSTREAMS,N1,Q) *   XNEG(I,1:NSTREAMS,N1)
                L_HOM    = DOT_PRODUCT ( LCON(1:NSTREAMS,N1), L_XPOS(I,1:NSTREAMS,N1,Q) ) &
                         + DOT_PRODUCT ( MCON(1:NSTREAMS,N1), HELP1 )
                COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
              ENDDO
            ENDDO

!  non-zero variations if LAYER_TO_VARY is an active layer above N
!    Get the linearized beam solution for the next layer
!  .. contributions from beam solutions on both sides.

          ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

            CALL LP_BEAMSOLUTION_NNEK &
             ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
               NSTREAMS, NSTREAMS_2, N1, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM, & ! Input, Numbers
               DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
               INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
               LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
               T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
               GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
               L_WUPPER, L_WLOWER )                                    ! Output

!mick eff 3/22/2017
            CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS_2
            DO Q = 1, N_LAYER_WFS
              COLTEL2_WF(CNS1:CNS2,Q) = L_WUPPER(1:NSTREAMS_2,N1,Q) - L_WLOWER(1:NSTREAMS_2,N,Q)
            ENDDO

!  Finish different types of boundary condition linearizations

          ENDIF

!  End loop over intermediate active layer boundaries

        ENDDO

!  Final boundary, bottom of lowest active layer
!  ---------------------------------------------

        NS = NLAYERS_TEL
        N  = ACTIVE_LAYERS(NS)      
        C0 = (NS-1)*NSTREAMS_2 + NSTREAMS

!  Get the linearized downward-reflected term. New 5/24/16 Rob Fix
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, each Fourier

        CALL LP_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, FOURIER,                 & ! Input
            LAYER_TO_VARY, N_LAYER_WFS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,  & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  If varying layer is same as last active telescoped layer, need full solutions

        IF ( N .EQ. LAYER_TO_VARY ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                CM = C0 + I
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
                HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                       - R2_L_HOMP(I,1:NSTREAMS,Q)
                HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
                L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
                COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
              ENDDO
            ENDDO
          ELSE
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                CM = C0 + I
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q)
                HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N)
                HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q)
                L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
                COLTEL2_WF(CM,Q) = - L_BEAM - L_HOM
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is below last active telescoped layer, need only the reflected solutions

        ELSE IF ( N .LT. LAYER_TO_VARY ) THEN
          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                CM = C0 + I
                I1 = I + NSTREAMS
                L_BEAM = R2_L_BEAM(I,Q)
                L_HOM  = DOT_PRODUCT(LCON(1:NSTREAMS,N),R2_L_HOMP(I,1:NSTREAMS,Q)) + &
                         DOT_PRODUCT(MCON(1:NSTREAMS,N),R2_L_HOMM(I,1:NSTREAMS,Q))
                COLTEL2_WF(CM,Q) = L_BEAM + L_HOM
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is above last telescoped layer, need reflected solutions
!     and beam solution linearizations propagated downwards

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                CM = C0 + I
                L_BEAM = - L_WLOWER(I1,N,Q) + R2_L_BEAM(I,Q)
                L_HOM  = DOT_PRODUCT(LCON(1:NSTREAMS,N),R2_L_HOMP(I,1:NSTREAMS,Q)) + &
                         DOT_PRODUCT(MCON(1:NSTREAMS,N),R2_L_HOMM(I,1:NSTREAMS,Q))
                COLTEL2_WF(CM,Q) = L_BEAM + L_HOM 
              ENDDO
            ENDDO
          ELSE
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                CM = C0 + I
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q)
                COLTEL2_WF(CM,Q) = - L_BEAM
              ENDDO
            ENDDO
          ENDIF

        ENDIF

 !  Add direct beam solution. This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              BEAM = DIRECT_BEAM(I,IBEAM)
              L_BEAM = - BEAM * ( L_DELTAU_VERT(Q,LAYER_TO_VARY) * DELTAU_SLANT(NLAYERS,LAYER_TO_VARY,IBEAM) )
              COLTEL2_WF(CM,Q) = COLTEL2_WF(CM,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
            ENDDO
          ENDDO
        ENDIF

!  Continuation point for the single-active-layer case

      ELSE

!  zero column vector

        SCOL2_WF = ZERO

!  top of active layer
!  -------------------

        NS = 1
        N = ACTIVE_LAYERS(NS)

!  If single active layer = layer that is varying,
!       then require homogeneous and beam solution linearizations

        IF ( LAYER_TO_VARY .EQ. N ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              L_BEAM = - L_WUPPER(I,N,Q)
              HELP1  = L_XPOS(I,1:NSTREAMS,N,Q)
              HELP2  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XNEG(I,1:NSTREAMS,N,Q) &
                     + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XNEG(I,1:NSTREAMS,N)
              L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                     + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
              SCOL2_WF(I,Q) = L_BEAM - L_HOM
            ENDDO
          ENDDO

!  Otherwise if varying layer is above active layer, there are beam
!  solution contributions propagated downwards - find these by calling
!  the appropriate solution module = L_BEAMSOLUTION_NNEK

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          CALL LP_BEAMSOLUTION_NNEK &
           ( DO_LAYER_SCATTERING, DO_PLANE_PARALLEL, TAYLOR_ORDER,           & ! Input, Flags and order
             NSTREAMS, NSTREAMS_2, N, LAYER_TO_VARY, N_LAYER_WFS, M, IBEAM,  & ! Input, Numbers
             DELTAU_VERT, SOLARBEAM_CUTOFF,                        & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
             LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
             T_DELT_EIGEN, XPOS, GFUNC_UP, GFUNC_DN,               & ! Input, Homogeneous/Green
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                    ! Output

!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(1:NSTREAMS,Q) = - L_WUPPER(1:NSTREAMS,N,Q)
          ENDDO

        ENDIF

!  Bottom of active layer
!  ----------------------

!  If active layer is varying layer, need full calculation.
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  get the linearized downward-reflected term. New 5/24/16 Rob Fix
!    -- 2/28/21. Version 3.8.3. BRDF_F array defined locally, each Fourier

        CALL LP_BVPTEL_SURFACE_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, N, FOURIER,                 & ! Input
            LAYER_TO_VARY, N_LAYER_WFS, SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F,  & ! Input
            T_DELT_EIGEN, L_T_DELT_EIGEN, T_DELT_DISORDS, L_T_DELT_DISORDS,    & ! Input
            XPOS, XNEG, WLOWER, L_XPOS, L_XNEG, L_WLOWER,                      & ! Input
            R2_L_BEAM, R2_L_HOMP, R2_L_HOMM, CUMTRANS, L_CUMTRANS )              ! Output 

!  If the single telescoped layer is the one that is varying

        IF ( N .EQ. LAYER_TO_VARY ) THEN
          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
                HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N) &
                       - R2_L_HOMP(I,1:NSTREAMS,Q)
                HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q) - R2_L_HOMM(I,1:NSTREAMS,Q)
                L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
                SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
              ENDDO
            ENDDO
          ELSE
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q)
                HELP1  =   T_DELT_EIGEN(1:NSTREAMS,N)   * L_XPOS(I1,1:NSTREAMS,N,Q) &
                       + L_T_DELT_EIGEN(1:NSTREAMS,N,Q) *   XPOS(I1,1:NSTREAMS,N)
                HELP2  = L_XNEG(I1,1:NSTREAMS,N,Q)
                L_HOM  = DOT_PRODUCT ( LCON(1:NSTREAMS,N), HELP1 ) &
                       + DOT_PRODUCT ( MCON(1:NSTREAMS,N), HELP2 )
                SCOL2_WF(I1,Q) = - L_BEAM - L_HOM
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is below single telescoped layer, need only the reflected solutions

        ELSE IF ( N .LT. LAYER_TO_VARY ) THEN
          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                L_BEAM = R2_L_BEAM(I,Q)
                L_HOM  = DOT_PRODUCT(LCON(1:NSTREAMS,N),R2_L_HOMP(I,1:NSTREAMS,Q)) + &
                         DOT_PRODUCT(MCON(1:NSTREAMS,N),R2_L_HOMM(I,1:NSTREAMS,Q))
                SCOL2_WF(I1,Q) = L_BEAM + L_HOM
              ENDDO
            ENDDO
          ENDIF

!   If varying layer is above single telescoped layer, need reflected solutions
!     and beam solution linearizations propagated downwards

        ELSE IF ( LAYER_TO_VARY .LT. N ) THEN

          IF ( DO_INCLUDE_SURFACE ) THEN
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                L_BEAM = - L_WLOWER(I1,N,Q) + R2_L_BEAM(I,Q)
                L_HOM  = DOT_PRODUCT(LCON(1:NSTREAMS,N),R2_L_HOMP(I,1:NSTREAMS,Q)) + &
                         DOT_PRODUCT(MCON(1:NSTREAMS,N),R2_L_HOMM(I,1:NSTREAMS,Q))
                SCOL2_WF(I1,Q) = L_BEAM + L_HOM 
              ENDDO
            ENDDO
          ELSE
!mick eff 3/22/2017
            DO Q = 1, N_LAYER_WFS
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                L_BEAM = L_WLOWER(I1,N,Q)
                SCOL2_WF(I1,Q) = - L_BEAM
              ENDDO
            ENDDO
          ENDIF

        ENDIF

 !  Add direct beam solution. This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

        IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
!mick eff 3/22/2017
          DO Q = 1, N_LAYER_WFS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              BEAM = DIRECT_BEAM(I,IBEAM) 
              L_BEAM = - BEAM * ( L_DELTAU_VERT(Q,LAYER_TO_VARY) * DELTAU_SLANT(NLAYERS,LAYER_TO_VARY,IBEAM) )
              SCOL2_WF(I1,Q) = SCOL2_WF(I1,Q) + L_BEAM * CUMTRANS(I) + L_CUMTRANS(I,Q) * BEAM
            ENDDO
          ENDDO
        ENDIF

!  End layer clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE LP_BVPTEL_COLUMN_SETUP

!  End

end module lidort_lp_bvproblem_m
