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

!    ##########################################################################
!    #                                                                        #
!    #  Collection of LAPACK Software Routines used in LIDORT Version 3.8.1   #
!    #                                                                        #
!    #    LAPACK Software comes with the following Modified BSD License       #
!    #                                                                        #
!    # Copyright (c) 1992-2013 The University of Tennessee and The University #
!    #                         of Tennessee Research Foundation.  All rights  #
!    #                         reserved.                                      #
!    # Copyright (c) 2000-2013 The University of California Berkeley. All     #
!    #                         rights reserved.                               #
!    # Copyright (c) 2006-2013 The University of Colorado Denver.  All rights #
!    #                         reserved.                                      #
!    #                                                                        #
!    # $COPYRIGHT$                                                            #
!    #                                                                        #
!    # Additional copyrights may follow                                       #
!    #                                                                        #
!    # $HEADER$                                                               #
!    #                                                                        #
!    # Redistribution and use in source and binary forms, with or without     #
!    # modification, are permitted provided that the following conditions are #
!    # met:                                                                   #
!    #                                                                        #
!    # - Redistributions of source code must retain the above copyright       #
!    #   notice, this list of conditions and the following disclaimer.        #
!    #                                                                        #
!    # - Redistributions in binary form must reproduce the above copyright    #
!    #   notice, this list of conditions and the following disclaimer listed  #
!    #   in this license in the documentation and/or other materials          #
!    #   provided with the distribution.                                      #
!    #                                                                        #
!    # - Neither the name of the copyright holders nor the names of its       #
!    #   contributors may be used to endorse or promote products derived from #
!    #   this software without specific prior written permission.             #
!    #                                                                        #
!    # The copyright holders provide no reassurances that the source code     #
!    # provided does not infringe any patent, copyright, or any other         #
!    # intellectual property rights of third parties.  The copyright holders  #
!    # disclaim any liability to any recipient for claims brought against     #
!    # recipient by any third party for infringement of that parties          #
!    # intellectual property rights.                                          #
!    #                                                                        #
!    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    #
!    # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      #
!    # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  #
!    # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT   #
!    # OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  #
!    # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       #
!    # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  #
!    # DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  #
!    # THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    #
!    # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  #
!    # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   #
!    #                                                                        #
!    ##########################################################################

! ###############################################################
! #                                                             #
! #  This is LIDORT_AUX.f90. This is a collection of utility    #
! #  modules compiled from a variety of sources. This is the    #
! #  only part of the LIDORT package not written by the Author. #
! #  The subroutines in LIDORT_AUX are listed below with their  #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      lcstring:              E. Mikusch, DLR 1994            #
! #      len_string:            E. Mikusch, DLR 1994            #
! #      lidort_read_error:     M. Christi, 2017                #
! #      lidort_write_status:   R.J.D. Spurr, RTS, 2005         #
! #      asymtx:                adapted by RJDS 1999 from the   #
! #                             DISORT module ASYMTX            #
! #      getquad2               M. Christi, 2017                #
! #      issort                 M. Christi, 2017                #
! #      rssort                 M. Christi, 2017                #
! #      rssort_idx             M. Christi, 2017                #
! #      cfplgarr:              T.P. Kurosu, 1997 (adapted      #
! #                             Numerical Recipes).             #
! #      xfindpar:            J. Lavagnino, SAO, 1991           #
! #      gfindpar:            J. Lavagnino, SAO, 1991           #
! #      FINDPAR_ERROR            R. Spurr, 2001                #
! #                                                             #
! #    LAPACK ROUTINES used in LIDORT                           #
! #       ( see separaate licensing statement)                  #
! #                                                             #
! #      dgbtrf:                  LAPACK                        #
! #      dgbtrs:                  LAPACK                        #
! #      dgetrf:                  LAPACK                        #
! #      dgetrs:                  LAPACK                        #
! #      dgbtf2:                  LAPACK                        #
! #      dlaswp:                  LAPACK                        #
! #      ilaenv:                  LAPACK                        #
! #      lsame:                   LAPACK                        #
! #      xerbla:                  LAPACK                        #
! #      dgetf2:                  LAPACK                        #
! #      dcopy:                   LAPACK                        #
! #      dgemm:                   LAPACK                        #
! #      dgemv:                   LAPACK                        #
! #      dger:                    LAPACK                        #
! #      dscal:                   LAPACK                        #
! #      dswap:                   LAPACK                        #
! #      dtbsv:                   LAPACK                        #
! #      dtrsm:                   LAPACK                        #
! #      idamax:                  LAPACK                        #
! #                                                             #
! ###############################################################

!  2/28/21. Version 3.8.3. Changed the call to ASYMTX to include tolerance variable

module lidort_aux_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk,ZERO,QUARTER,HALF,ONE,TWO,PIE

!  Only the Top-level LAPACK routines are public

   PRIVATE
   PUBLIC ::    LEN_STRING,                      & ! Bookkeeping
                LIDORT_READ_ERROR,               & ! Error Handling
                LIDORT_WRITE_STATUS,             & ! Error Handling
                ASYMTX, GETQUAD2, ISSORT,RSSORT, & ! Numerical (Not LAPACK)
                RSSORT_IDX, CFPLGARR,            & ! Numerical (Not LAPACK)
                DGBTRF, DGBTRS, DGETRF, DGETRS,  & ! Lapack Top-level
                GFINDPAR, FINDPAR_ERROR            ! Findpar

contains

SUBROUTINE LCSTRING ( STRING )

!  Convert a string to lower-case.  (ASCII character set assumed.)

      implicit none

!  Modified argument.

      CHARACTER (LEN=*), intent(inout)   :: STRING

!  Local variables.

      INTEGER         :: CCC
      INTEGER         :: III
      INTEGER         :: OFFSET

!***********************************************************************

      OFFSET = ICHAR ( 'a' ) - ICHAR ( 'A' )
      DO 10 III = 1, LEN ( STRING )
        CCC = ICHAR ( STRING ( III : III ) )
        IF ( CCC .GE. ICHAR ( 'A' ) .AND.  CCC .LE. ICHAR ( 'Z' ) ) THEN
          STRING ( III : III ) = CHAR ( CCC + OFFSET )
        END IF
10    CONTINUE

END SUBROUTINE LCSTRING

!

INTEGER FUNCTION LEN_STRING(STR)

      implicit none

!     ----------------------------------------------------------------
!       Length of string.
!       This function determines the length of the string contained in 
!       the character variable. Therefore it searches the last non-blank
!       character and returnes the position of this value as the length.
!       If a character variable contains blanks only, 0 is returned.
!        2 Feb 1994, Eberhard Mikusch 
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!     input argument
!     ----------------------------------------------------------------
      CHARACTER (LEN=*), INTENT(IN) :: STR

!     ----------------------------------------------------------------
!     get last non-blank character
!     ----------------------------------------------------------------
      LEN_STRING = LEN(STR)

      DO WHILE ((LEN_STRING .GT.0) .AND.  &
               (STR(LEN_STRING:LEN_STRING) .EQ. ' '))
        LEN_STRING = LEN_STRING - 1
      END DO 

!     ----------------------------------------------------------------
!     return
!     ----------------------------------------------------------------
      RETURN
END FUNCTION LEN_STRING

!

SUBROUTINE LIDORT_READ_ERROR ( ERRORFILE, LIDORT_InputStatus )

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : LIDORT_ERRUNIT
      USE LIDORT_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  LIDORT file inputs status

      TYPE(LIDORT_Input_Exception_Handling) :: LIDORT_InputStatus

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = LIDORT_ERRUNIT

!  Write LIDORT configuration file read errors to LIDORT error file 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL:   Wrong input from LIDORT input file-read'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',LIDORT_InputStatus%TS_NINPUTMESSAGES
      DO N = 1, LIDORT_InputStatus%TS_NINPUTMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',&
          ADJUSTL(TRIM(LIDORT_InputStatus%TS_INPUTMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',&
          ADJUSTL(TRIM(LIDORT_InputStatus%TS_INPUTACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Read-input fail: Look at file ' // TRIM(ERRORFILE)
      STOP

END SUBROUTINE LIDORT_READ_ERROR

!

SUBROUTINE LIDORT_WRITE_STATUS &
        ( ERRORFILE, ERRORUNIT, OPENFILEFLAG, LIDORT_Status )

!mick mod 3/22/2017 - added variable LIDORT_ABORT
!                   - enhanced internal error handling for fatal cases

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_MESSAGES,   LIDORT_SUCCESS, &
                                LIDORT_WARNING, LIDORT_SERIOUS
      USE LIDORT_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file and unit (will only opened if overall status not successful)

      CHARACTER (LEN=*), intent(in) :: ERRORFILE
      INTEGER, intent(in)           :: ERRORUNIT

      LOGICAL, intent(inout)        :: OPENFILEFLAG

!  LIDORT status

      TYPE(LIDORT_Exception_Handling), intent(in) :: LIDORT_Status

!  Local variables

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Help

      INTEGER          :: N, W
      LOGICAL          :: STATUS_OVERALL, LIDORT_ABORT

!  Copy inputs to local variables

      STATUS_INPUTCHECK               = LIDORT_Status%TS_STATUS_INPUTCHECK
      NCHECKMESSAGES                  = LIDORT_Status%TS_NCHECKMESSAGES
      CHECKMESSAGES(0:NCHECKMESSAGES) = LIDORT_Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES)
      ACTIONS(0:NCHECKMESSAGES)       = LIDORT_Status%TS_ACTIONS(0:NCHECKMESSAGES)

!mick fix 3/22/2017 - added IF condition and ELSE section
      IF ( STATUS_INPUTCHECK .NE. LIDORT_SERIOUS ) THEN
        STATUS_CALCULATION = LIDORT_Status%TS_STATUS_CALCULATION
        MESSAGE            = LIDORT_Status%TS_MESSAGE
        TRACE_1            = LIDORT_Status%TS_TRACE_1
        TRACE_2            = LIDORT_Status%TS_TRACE_2
        TRACE_3            = LIDORT_Status%TS_TRACE_3

        STATUS_OVERALL = ( STATUS_INPUTCHECK  .NE. LIDORT_SUCCESS ) .or. &
                         ( STATUS_CALCULATION .NE. LIDORT_SUCCESS )
      ELSE
        STATUS_OVERALL = .true.
      ENDIF

!  Open if flagged

      IF ( .not.OPENFILEFLAG .and. STATUS_OVERALL ) THEN
        OPENFILEFLAG = .true.
        OPEN(ERRORUNIT, FILE = trim(ERRORFILE), STATUS = 'unknown' )
      ENDIF

!  Define some local variables

      W = ERRORUNIT
      LIDORT_ABORT = .false.

!  Check status and write messages

      IF ( STATUS_INPUTCHECK .EQ. LIDORT_WARNING ) THEN
        WRITE(W,*)' WARNING: Suspicious input from LIDORT input check'
        WRITE(W,*)'  ------ LIDORT has executed with internal defaults'
        WRITE(W,*)'  ------ Here are the messages and actions '
        WRITE(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',adjustl(trim(CHECKMESSAGES(N)))
          WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',adjustl(trim(ACTIONS(N)))
        ENDDO
      ELSE IF ( STATUS_INPUTCHECK .EQ. LIDORT_SERIOUS ) THEN
        WRITE(W,*)' FATAL:   Wrong input from LIDORT input check '
        WRITE(W,*)'  ------ LIDORT will not execute'
        WRITE(W,*)'  ------ Here are the messages and actions '
        WRITE(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',adjustl(trim(CHECKMESSAGES(N)))
          WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',adjustl(trim(ACTIONS(N)))
        ENDDO

        !Set abort flag and alert the user at the screen
        LIDORT_ABORT = .true.
        WRITE(*,*)
        WRITE(*,*)'LIDORT input abort'
      ENDIF

!  Model calculation status

      IF ( STATUS_INPUTCHECK .NE. LIDORT_SERIOUS ) THEN
        IF ( STATUS_CALCULATION .NE. LIDORT_SUCCESS ) THEN
          WRITE(W,*)' FATAL: LIDORT execution failure'
          WRITE(W,*)'  -----Here is the message and traces : '
          WRITE(W,'(A,A)')' Message : ', adjustl(trim(MESSAGE))
          WRITE(W,'(A,A)')' Trace 1 : ', adjustl(trim(TRACE_1))
          WRITE(W,'(A,A)')' Trace 2 : ', adjustl(trim(TRACE_2))
          WRITE(W,'(A,A)')' Trace 3 : ', adjustl(trim(TRACE_3))

          !Set abort flag and alert the user at the screen
          LIDORT_ABORT = .true.
          WRITE(*,*)
          WRITE(*,*)'LIDORT calculation abort'
        ENDIF
      ENDIF

!  If fatal error, abort

      IF ( LIDORT_ABORT ) THEN
        CLOSE(ERRORUNIT)
        WRITE(*,*)
        WRITE(*,'(/1x,a)')'Main program failed to execute, ' // &
          'error messages in ' // trim(ERRORFILE)
        STOP
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_WRITE_STATUS

!

SUBROUTINE ASYMTX ( TOL, AAD, M, IA, IEVEC, EVECD, EVALD, IER, WKD, &
                    MESSAGE, BAD_STATUS )

!  2/28/21. Version 3.8.3. Add the tolerance variable

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, ZERO, ONE

      implicit none

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.

!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.

!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the parlett-reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.

!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

!   I N P U T    V A R I A B L E S:

!        AAD  :  input asymmetric matrix, destroyed after solved
!        M    :  order of  A
!       IA    :  first dimension of  A
!    IEVEC    :  first dimension of  EVECD

!   O U T P U T    V A R I A B L E S:

!       EVECD :  (unnormalized) eigenvectors of  A 
!                   ( column J corresponds to EVALD(J) )

!       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

      LOGICAL      , intent(out) :: BAD_STATUS
      CHARACTER*(*), intent(out) :: MESSAGE

!   S C R A T C H   V A R I A B L E S:

!       WKD    :  WORK AREA ( DIMENSION AT LEAST 2*M )

!  input/output arguments
!    - 2/28/21. Version 3.8.3. Add the tolerance variable

      INTEGER, intent(in)   :: M, IA, IEVEC
      REAL(fpk), intent(in) :: TOL
      INTEGER, intent(out)  :: IER

      REAL(fpk), intent(inout) :: AAD(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(inout) :: WKD(MAXSTREAMS_2)

      REAL(fpk), intent(out)   :: EVALD(MAXSTREAMS)
      REAL(fpk), intent(out)   :: EVECD(MAXSTREAMS,MAXSTREAMS)

!  local variables (explicit declaration)

      LOGICAL          :: NOCONV, NOTLAS
      INTEGER          :: I, J, L, K, KKK, LLL
      INTEGER          :: N, N1, N2, IN, LB, KA, II
      !REAL(fpk)        :: C1, C2, C3, C4, C5, C6
      !DATA               C1 / 0.4375D0 /
      !DATA               C2 / 0.5D0 /
      !DATA               C3 / 0.75D0 /
      !DATA               C4 / 0.95D0 /
      !DATA               C5 / 16.0D0 /
      !DATA               C6 / 256.0D0 /
      REAL(fpk), PARAMETER :: C1=0.4375D0, &
                              C2=0.5D0,    &
                              C3=0.75D0,   &
                              C4=0.95D0,   &
                              C5=16.0D0,   &
                              C6=256.0D0
      REAL(fpk)        :: DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      REAL(fpk)        :: REPL, COL, ROW, SCALE, T, X, Z, S, Y, UU, VV

!  these variables by R. Spurr, Removes GOTO statements
!     Introduced 08 Feburary 2010, RT SOLUTIONS Inc.

      LOGICAL          :: B70, B120, B200, B380, B390, B400, B460
      LOGICAL          :: TBLOCK, HESSENBERG_FORM, B490, B520

!  output status

      IER        = 0
      BAD_STATUS = .FALSE.
      MESSAGE    = ' '

!       Here change to bypass D1MACH:
!  2/28/21. Version 3.8.3. Removed this section
!      TOL = 0.0000001
!       TOL = 1.0D-20
!        TOL = D1MACH(4)

      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) THEN
        MESSAGE = 'ASYMTX--bad input variable(s)'
        BAD_STATUS = .TRUE.
        RETURN
      ENDIF
!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES

      IF ( M.EQ.1 )  THEN
         EVALD(1) = AAD(1,1)
         EVECD(1,1) = ONE
         RETURN
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AAD(1,1) - AAD(2,2) )**2 + 4.0D0*AAD(1,2)*AAD(2,1)
         IF ( DISCRI.LT.ZERO ) THEN
           MESSAGE = 'ASYMTX--COMPLEX EVALS IN 2X2 CASE'
           BAD_STATUS = .TRUE.
           RETURN
         ENDIF
         SGN = ONE
         IF ( AAD(1,1).LT.AAD(2,2) )  SGN = - ONE
         EVALD(1) = 0.5D0*( AAD(1,1) + AAD(2,2) + SGN*DSQRT(DISCRI) )
         EVALD(2) = 0.5D0*( AAD(1,1) + AAD(2,2) - SGN*DSQRT(DISCRI) )
         EVECD(1,1) = ONE
         EVECD(2,2) = ONE
         IF ( AAD(1,1).EQ.AAD(2,2) .AND. &
              (AAD(2,1).EQ.ZERO.OR.AAD(1,2).EQ.ZERO) ) THEN
            RNORM = DABS(AAD(1,1))+DABS(AAD(1,2))+ &
                     DABS(AAD(2,1))+DABS(AAD(2,2))
            W = TOL * RNORM
            EVECD(2,1) = AAD(2,1) / W
            EVECD(1,2) = - AAD(1,2) / W
         ELSE
            EVECD(2,1) = AAD(2,1) / ( EVALD(1) - AAD(2,2) )
            EVECD(1,2) = AAD(1,2) / ( EVALD(2) - AAD(1,1) )
         ENDIF
         RETURN
      END IF

!                                        ** INITIALIZE OUTPUT VARIABLES
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE


!                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
!                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
!                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M

      TBLOCK = .TRUE.    
      DO WHILE (TBLOCK)

         KKK = K 
         B70 = .TRUE.
         J = KKK + 1

         DO WHILE (B70.and.TBLOCK)

            J = J - 1
            TBLOCK = (J.GT.1)

            ROW = ZERO
            DO 40 I = 1, K
               IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )
40          CONTINUE

            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               B70 = .false.
  
            ENDIF

         ENDDO
      ENDDO

!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT

      TBLOCK = .TRUE.    
      DO WHILE (TBLOCK)

         LLL = L
         B120 = .TRUE.
         J = LLL - 1

         DO WHILE (B120.and.TBLOCK)

            J = J + 1
            TBLOCK = (J.LT.K)

            COL = ZERO

            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )
90          CONTINUE

            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               B120 = .false.
             ENDIF

           ENDDO

       ENDDO

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!       pause'after 120'

!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K

      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE

      B200 = .true.

      DO WHILE ( B200 )

         noconv = .false.

         DO 200 I = L, K

            COL = ZERO
            ROW = ZERO

            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + DABS( AAD(J,I) )
                  ROW = ROW + DABS( AAD(I,J) )
               END IF
150         CONTINUE

            F = ONE
            G = ROW / C5
            H = COL + ROW

            DO WHILE ( COL .LT. G )
               F   = F * C5
               COL = COL * C6
            ENDDO

            G = ROW * C5
            DO WHILE ( COL .GE. G )
               F   = F / C5
               COL = COL / C6
            ENDDO

            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
           
200      CONTINUE

         B200 = ( noconv.eqv..true.)

      ENDDO

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!       pause'before 350'

!  ** IS -A- ALREADY IN HESSENBERG FORM?

      HESSENBERG_FORM = ( K-1 .LT. L+1 )

      IF ( .not. HESSENBERG_FORM ) THEN

!   ** TRANSFER -A- TO A HESSENBERG FORM

         DO 290 N = L+1, K-1
            H        = ZERO
            WKD(N+M) = ZERO
            SCALE    = ZERO
!                                                        ** SCALE COLUMN
            DO 210 I = N, K
               SCALE = SCALE + DABS(AAD(I,N-1))
210         CONTINUE
            IF ( SCALE.NE.ZERO ) THEN
               DO 220 I = K, N, -1
                  WKD(I+M) = AAD(I,N-1) / SCALE
                  H = H + WKD(I+M)**2
220            CONTINUE
               G = - SIGN( DSQRT(H), WKD(N+M) )
               H = H - WKD(N+M) * G
               WKD(N+M) = WKD(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
               DO 250 J = N, M
                  F = ZERO
                  DO 230  I = K, N, -1
                     F = F + WKD(I+M) * AAD(I,J)
230               CONTINUE
                  DO 240 I = N, K
                     AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240               CONTINUE
250            CONTINUE
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
               DO 280 I = 1, K
                  F = ZERO
                  DO 260  J = K, N, -1
                     F = F + WKD(J+M) * AAD(I,J)
260               CONTINUE
                  DO 270 J = N, K
                     AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270               CONTINUE
280            CONTINUE
               WKD(N+M)  = SCALE * WKD(N+M)
               AAD(N,N-1) = SCALE * G
            END IF
290      CONTINUE

         DO 340  N = K-2, L, -1
            N1 = N + 1
            N2 = N + 2
            F  = AAD(N1,N)
            IF ( F.NE.ZERO ) THEN
               F  = F * WKD(N1+M)
               DO 300 I = N2, K
                  WKD(I+M) = AAD(I,N)
300            CONTINUE
               IF ( N1.LE.K ) THEN
                  DO 330 J = 1, M
                     G = ZERO
                     DO 310 I = N1, K
                        G = G + WKD(I+M) * EVECD(I,J)
310                  CONTINUE
                     G = G / F
                     DO 320 I = N1, K
                        EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320                  CONTINUE
330               CONTINUE
               END IF
            END IF
340      CONTINUE

!  End clause for conversion to Hessenberg Form

      ENDIF

!      do i = 1, 4
!       write(*,'(4f10.6)')(AAD(I,J),J=1,4)
!      enddo
!      do i = 1, 4
!       write(*,'(4f10.6)')(EVECD(I,J),J=1,4)
!      enddo
!       pause'before 350'

!  Set first Eigenvalues

      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + DABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)
370   CONTINUE
      N = K
      T = ZERO

! #################################################################
! #################################################################
! ################################################################

!    ** SEARCH FOR NEXT EIGENVALUES

      B380 = ( N.GE.L )
      DO WHILE (B380 )

         IN = 0
         N1 = N - 1
         N2 = N - 2

!  ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
 
         B390 = .true.
         DO WHILE ( B390 )

! 400  Block

            I = L - 1
            B400 = .true.
            do while ( B400 .and. I.lt.N )
               I = I + 1
               LB = N+L - I
               IF ( LB.EQ.L ) B400 = .false.
               if ( B400 ) then
                 S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
                 IF ( S.EQ.ZERO ) S = RNORM
                 IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) B400 = .false.
               endif
            enddo

!   ** ONE EIGENVALUE FOUND

            X = AAD(N,N)      
            IF ( LB.EQ.N ) THEN
               AAD(N,N)  = X + T
               EVALD(N) = AAD(N,N)
               N = N1
               B390 = .false.
               B380 = ( N.GE.L )
            END IF

!  Other wise ** TWO EIGENVALUES FOUND

            IF ( B390 ) THEN

               Y = AAD(N1,N1)
               W = AAD(N,N1) * AAD(N1,N)

              IF ( LB.EQ.N1 ) THEN
                  P = (Y-X) * C2
                  Q = P**2 + W
                  Z = DSQRT( DABS(Q) )
                  AAD(N,N) = X + T
                  X = AAD(N,N)
                  AAD(N1,N1) = Y + T
!                                        ** REAL PAIR
                  Z = P + SIGN(Z,P)
                  EVALD(N1) = X + Z
                  EVALD(N)  = EVALD(N1)
                  IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z
                  X = AAD(N,N1)
!                                  ** EMPLOY SCALE FACTOR IN CASE
!                                  ** X AND Z ARE VERY SMALL
                  R = SQRT( X*X + Z*Z )
                  P = X / R
                  Q = Z / R
!                                             ** ROW MODIFICATION
                  DO 420 J = N1, M
                     Z = AAD(N1,J)
                     AAD(N1,J) = Q * Z + P * AAD(N,J)
                     AAD(N,J)  = Q * AAD(N,J) - P * Z
420               CONTINUE
!                                             ** COLUMN MODIFICATION
                  DO 430 I = 1, N
                     Z = AAD(I,N1)
                     AAD(I,N1) = Q * Z + P * AAD(I,N)
                     AAD(I,N)  = Q * AAD(I,N) - P * Z
430               CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
                  DO 440 I = L, K
                     Z = EVECD(I,N1)
                     EVECD(I,N1) = Q * Z + P * EVECD(I,N)
                     EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440               CONTINUE

                  N = N2
                  B390 = .false.
                  B380 = ( N.GE.L )
               END IF

            ENDIF

!  Carry on

            if ( B390 ) THEN

!    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR & return
!    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE

               IF ( IN.EQ.30 ) THEN
                  IER = N
                  RETURN
               END IF

!     ** FORM SHIFT

               IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
                   T = T + X
                   DO 450 I = L, N
                      AAD(I,I) = AAD(I,I) - X
450                CONTINUE
                   S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
                   X = C3 * S
                   Y = X
                   W = - C1 * S**2
                END IF

                IN = IN + 1

!  ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

                B460 = .true.
                J = LB - 1 
                DO WHILE (B460.and.J.LT.N2)
                   J = J + 1
                   I = N2+LB - J
                   Z = AAD(I,I)
                   R = X - Z
                   S = Y - Z
                   P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
                   Q = AAD(I+1,I+1) - Z - R - S
                   R = AAD(I+2,I+1)
                   S = DABS(P) + DABS(Q) + DABS(R)
                   P = P / S
                   Q = Q / S
                   R = R / S
                   IF ( I.EQ.LB ) B460 = .false.
                   if ( B460 ) then
                      UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
                      VV = DABS(P) * &
                       (DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))
                     IF ( UU .LE. TOL*VV ) B460 = .false.
                   endif
                ENDDO

!  B470 was here

               AAD(I+2,I) = ZERO
               DO 480 J = I+3, N
                  AAD(J,J-2) = ZERO
                  AAD(J,J-3) = ZERO
480            CONTINUE

!   ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

               KA = I - 1
               B520 = .true.

               DO WHILE (B520.and.KA.lt.N1)

                  KA = KA + 1 
                  NOTLAS = KA.NE.N1
                  B490 = .true.

                  IF ( KA.EQ.I ) THEN
                     S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
                     IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
                  ELSE
                     P = AAD(KA,KA-1)
                     Q = AAD(KA+1,KA-1)
                     R = ZERO
                     IF ( NOTLAS ) R = AAD(KA+2,KA-1)
                     X = DABS(P) + DABS(Q) + DABS(R)
                     IF ( X.EQ.ZERO ) THEN
                        B490 = .false.
                     ELSE
                        P = P / X
                        Q = Q / X
                        R = R / X
                        S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
                        AAD(KA,KA-1) = - S * X
                     END IF
                  ENDIF

!  Only do reminder if set

                  IF ( B490 ) then

                     P = P + S
                     X = P / S
                     Y = Q / S
                     Z = R / S
                     Q = Q / P
                     R = R / P
!                                                    ** ROW MODIFICATION
                     DO 490 J = KA, M
                        P = AAD(KA,J) + Q * AAD(KA+1,J)
                        IF ( NOTLAS ) THEN
                           P = P + R * AAD(KA+2,J)
                           AAD(KA+2,J) = AAD(KA+2,J) - P * Z
                        END IF
                        AAD(KA+1,J) = AAD(KA+1,J) - P * Y
                        AAD(KA,J)   = AAD(KA,J)   - P * X
490                  CONTINUE
!                                                 ** COLUMN MODIFICATION
                     DO 500 II = 1, MIN0(N,KA+3)
                        P = X * AAD(II,KA) + Y * AAD(II,KA+1)
                        IF ( NOTLAS ) THEN
                           P = P + Z * AAD(II,KA+2)
                           AAD(II,KA+2) = AAD(II,KA+2) - P * R
                        END IF
                        AAD(II,KA+1) = AAD(II,KA+1) - P * Q
                        AAD(II,KA)   = AAD(II,KA) - P
500                  CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
                     DO 510 II = L, K
                        P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
                        IF ( NOTLAS ) THEN
                           P = P + Z * EVECD(II,KA+2)
                           EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
                        END IF
                        EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
                        EVECD(II,KA)   = EVECD(II,KA) - P
510                  CONTINUE

!  B490 clause

                  ENDIF

!  End B520 do while

               ENDDO

!  Clause for B390

            ENDIF

!  Finish loop for B390

         ENDDO

!  Finish Loop for B380

      ENDDO

!      do i = 1, 4
!       write(*,'(5F14.7)')(AAD(I,J),J=1,4), evald(i)
!      enddo
!      pause'after evals'

!  ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR

      IF ( RNORM.NE.ZERO ) THEN

         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE


!                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE

!                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 610  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600            CONTINUE
610         CONTINUE
         END IF

      END IF


!  Set egienvectors

      DO 625 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620      CONTINUE
625   CONTINUE



!                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE

      DO 660 I = K+1, M
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE

!                    
!  670 CONTINUE    ! Removed

!  Finish

      RETURN
END SUBROUTINE ASYMTX

!

SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!  Limits of interval

      REAL(FPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(FPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(FPK) :: IR, MR, NR
      REAL(FPK) :: MIDPT, SFAC
      REAL(FPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(FPK), PARAMETER :: QEPS = 1.0D-13

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,FPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,FPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_FPK*NR**3) &
             - ONE/(384.0_FPK*NR**4)*(39.0_FPK - 28.0_FPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,FPK) ; LPM2 = LPM1 ; LPM1 = LP
               LP = ((TWO*IR - ONE)*X*LPM1 - (IR - ONE)*LPM2)/IR
            ENDDO
            DLP_DX = NR*(X*LP - LPM1)/(X**2 - ONE)
            XOLD = X ; X = XOLD - LP/DLP_DX
            IF (ABS(X-XOLD) <= QEPS) EXIT
         ENDDO

!  Shift and scale the current root (and its symmetric counterpart) from the interval (-1,1)
!  to the interval (a,b).  Define their related weights (e.g. see Abramowitz & Stegan (1972)).
!  Note:
!  If (1) N is even or (2) N is odd and M /= N2, then ROOTS(M) and ROOTS(NM1) are unique.
!  If N is odd and M = N2, then M = NM1 and ROOTS(M) = ROOTS(NM1) are one and the same root.

         !On interval lower half: (a,midpt)
         ROOTS(M)   = MIDPT - SFAC*X
         WGTS(M)    = (TWO*SFAC)/((ONE - X**2)*DLP_DX**2)

         !On interval upper half: (midpt,b)
         NM1 = N - M + 1
         ROOTS(NM1) = MIDPT + SFAC*X
         WGTS (NM1) = WGTS(M)

      ENDDO

END SUBROUTINE GETQUAD2

!

SUBROUTINE ISSORT(N,ILIST)

!  Sorts a 1-D array of integers in ascending order using a simple selection
!  ("interchange") sort if necessary.  Can handle cases where elements in the
!  array appear more than once.

      IMPLICIT NONE

!  Dimension

      INTEGER, INTENT(IN) :: N

!  List to sort

      INTEGER, INTENT(INOUT) :: ILIST(N)

!  Local variables

      INTEGER :: I, ITEMP, J, SMALL_INDEX
      LOGICAL :: IN_ORDER

!  Do prechecks

      !(1) Case where N = 1
      IF (N == 1) RETURN

      !(2) Case where list is already in order
      IN_ORDER = .TRUE.
      DO I = 1, N-1
        IF (ILIST(I) > ILIST(I+1)) THEN
          IN_ORDER = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (IN_ORDER) RETURN

!  Sort

      DO I = 1, N-1
        SMALL_INDEX = I
        DO J = I+1, N
          IF (ILIST(J) < ILIST(SMALL_INDEX)) SMALL_INDEX = J
        ENDDO
        IF (SMALL_INDEX /= I) THEN
          ITEMP              = ILIST(I)
          ILIST(I)           = ILIST(SMALL_INDEX)
          ILIST(SMALL_INDEX) = ITEMP
        ENDIF
      ENDDO

END SUBROUTINE ISSORT

!

SUBROUTINE RSSORT(N,RLIST)

!  Sorts a 1-D array of reals in ascending order using a simple selection
!  ("interchange") sort if necessary.  Can handle cases where elements in the
!  array appear more than once.

      IMPLICIT NONE

!  Dimension

      INTEGER, INTENT(IN) :: N

!  List to sort

      REAL(FPK), INTENT(INOUT) :: RLIST(N)

!  Local variables

      INTEGER   :: I, J, SMALL_INDEX
      REAL(FPK) :: RTEMP
      LOGICAL   :: IN_ORDER

!  Do prechecks

      !(1) Case where N = 1
      IF (N == 1) RETURN

      !(2) Case where list is already in order
      IN_ORDER = .TRUE.
      DO I = 1, N-1
        IF (RLIST(I) > RLIST(I+1)) THEN
          IN_ORDER = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (IN_ORDER) RETURN

!  Sort

      DO I = 1, N-1
        SMALL_INDEX = I
        DO J = I+1, N
          IF (RLIST(J) < RLIST(SMALL_INDEX)) SMALL_INDEX = J
        ENDDO
        IF (SMALL_INDEX /= I) THEN
          RTEMP              = RLIST(I)
          RLIST(I)           = RLIST(SMALL_INDEX)
          RLIST(SMALL_INDEX) = RTEMP
        ENDIF
      ENDDO

END SUBROUTINE RSSORT

!

SUBROUTINE RSSORT_IDX(N,RLIST,IDX)

!  Takes a 1-D input array of reals RLIST and returns a 1-D array of integers
!  IDX whose entries are the indices of RLIST in an order needed to put the
!  elements of RLIST in ascending order.  Uses a simple selection ("interchange")
!  sort if necessary.  Can handle cases where elements in RLIST appear more than
!  once.

      IMPLICIT NONE

!  Dimension

      INTEGER, INTENT(IN) :: N

!  List whose indices may require sorting

      REAL(FPK), INTENT(IN) :: RLIST(N)

!  List of sorted indices

      INTEGER, INTENT(OUT) :: IDX(N)

!  Local variables

      INTEGER   :: I, ITEMP, J, SMALL_INDEX, SMALL_J
      LOGICAL   :: IN_ORDER

!  Initialize output

      DO I = 1, N
        IDX(I) = I
      END DO

!  Do prechecks

      !(1) Case where N = 1
      IF (N == 1) RETURN

      !(2) Case where list is already in order
      IN_ORDER = .TRUE.
      DO I = 1, N-1
        IF (RLIST(I) > RLIST(I+1)) THEN
          IN_ORDER = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (IN_ORDER) RETURN

!  Sort

      DO I = 1, N-1
        SMALL_INDEX = IDX(I)
        DO J = I+1, N
          IF (RLIST(IDX(J)) < RLIST(SMALL_INDEX)) THEN
            SMALL_INDEX = IDX(J)
            SMALL_J     = J
          END IF
        ENDDO
        IF (SMALL_INDEX /= IDX(I)) THEN
          ITEMP        = IDX(I)
          IDX(I)       = IDX(SMALL_J)
          IDX(SMALL_J) = ITEMP
        ENDIF
      ENDDO

END SUBROUTINE RSSORT_IDX

!

SUBROUTINE CFPLGARR (MAXMOM, L, M, X, CFPLG)

!-------- This Routine Calculates  -----------

! ###########################################################
! #                                                         #
! #          [(l - m)! / (l + m)!]^{1/2} P(l,m)(x)          #
! #                                                         #
! ###########################################################

!     Based on the PLG routine from the numerical receipes,
!     this revised routine calculates all legendre polynomials
!     from Plg(m,m) to Plg(l,m).
!     The maximum order for the Plg is OMAX

      IMPLICIT NONE

!     ---------------
!     Input Variables
!     ---------------

      Integer  , intent(in)  :: maxmom
      Integer  , intent(in)  :: l, m       ! order of legendre polynomial
      Real(fpk), intent(in)  :: x          ! abscissa

!     ---------------
!     Output Variables
!     ---------------

      Real(fpk), intent(out) ::  cfplg (0:maxmom) ! array of Plgs

!     ---------------
!     Local Variables
!     ---------------

      Integer    :: i, ll
      Real(fpk)  :: fact, fll, fmm, fmmp1, somx2

!     --- Compute F(M,M)  ---

      fmm = 1.0d0

      If ( m .gt. 0) Then
         somx2 = Dsqrt( (1.0d0-x) * (1.0d0+x) )
         fact = 1.0d0
         Do i = 1, m
            fmm = -fmm*fact*somx2 / Dsqrt( DBLE( (2*i)*(2*i-1) ) )
            fact = fact + 2.0d0
         End Do
      End If

      cfplg(m) = fmm

!     --- Compute F(M+1,M)  ---

      If (l .gt. m) Then
         fmmp1 = x * Dsqrt( DBLE( 2*m+1 ) ) * fmm
         cfplg(m+1) = fmmp1

!     --- Compute F(L,M) ---

         If ( l .gt. (m+1) ) Then
         
            Do ll = m+2, l
               fll = (fmmp1 * x * DBLE(2*ll-1) - fmm * DBLE(ll+m-1) &
                   * Dsqrt( (DBLE( ll-m-1 ) ) / DBLE( ll+m-1) ) )  &
                   / Dsqrt( DBLE( (ll+m)*(ll-m) ) )
               fmm   = fmmp1
               fmmp1 = fll
               cfplg(ll) = fmmp1
            End Do

         End If

      End If

      Return
END SUBROUTINE CFPLGARR

! 

SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER          :: INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER          :: IPIV( * )
      REAL(FPK)        :: AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL(FPK)      :: array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      !    !    !    +    +    +       !    !    !   u14  u25  u36
!      !    !    +    +    +    +       !    !   u13  u24  u35  u46
!      !   a12  a23  a34  a45  a56      !   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   !      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   !    !      m31  m42  m53  m64   !    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)        ::ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER          :: NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER          ::  I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                           JU, K2, KM, KV, NB, NW
      REAL(FPK)      ::   TEMP
!     ..
!     .. Local Arrays ..
      REAL(FPK)      ::   WORK13( LDWORK, NBMAX ), &
                          WORK31( LDWORK, NBMAX )

!  No Externals 
!     ..
!     .. External Functions ..
!      INTEGER            :: IDAMAX, ILAENV
!      EXTERNAL           IDAMAX, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL, &
!                         DSWAP, DTRSM, XERBLA

!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                   WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                   AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), 1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )                                   &
                    CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                               AB( KV, JJ+1 ), LDAB-1,             &
                               AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 ) INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 ) &
                 CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                             WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                           IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                             JB, J2, ONE, AB( KV+1, J ), LDAB-1,      &
                             AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2, &
                                JB, -ONE, AB( KV+1+JB, J ), LDAB-1,     &
                                AB( KV+1-JB, J+JB ), LDAB-1, ONE,       &
                                AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2, &
                                JB, -ONE, WORK31, LDWORK,               &
                                AB( KV+1-JB, J+JB ), LDAB-1, ONE,       &
                                AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                             JB, J3, ONE, AB( KV+1, J ), LDAB-1,       &
                             WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3, &
                                JB, -ONE, AB( KV+1+JB, J ), LDAB-1,      &
                                WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),   &
                                LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3, &
                                JB, -ONE, WORK31, LDWORK, WORK13,        &
                                LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 ) &
                 CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                              AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of DGBTRF
!
END SUBROUTINE DGBTRF


SUBROUTINE DGBTRS ( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER*1      ::  TRANS
      INTEGER          ::  INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER        ::    IPIV( * )
      REAL(FPK)      ::   AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by DGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) REAL(FPK)      :: array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) REAL(FPK)      :: array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)   , parameter   ::  ONE = 1.0D+0 
!     ..
!     .. Local Scalars ..
      LOGICAL            :: LNOTI, NOTRAN
      INTEGER            :: I, J, KD, L, LM
!     ..
!     .. External Functions ..
!      LOGICAL            :: LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
         LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN
!
      KD = KU + KL + 1
      LNOTI = KL.GT.0
!
      IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) &
                 CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB ) 
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
                         LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
!
         DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
                        AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
!
      ELSE
!
!        Solve A'*X = B.
!
         DO 30 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
                        LDAB, B( 1, I ), 1 )
   30    CONTINUE
!
!        Solve L'*X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
                           LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) &
                  CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of DGBTRS
!
END SUBROUTINE DGBTRS

!

SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER        ::    INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER        ::    IPIV( * )
      REAL(FPK)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(FPK)      :: array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)      ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            :: I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
!      INTEGER            :: ILAENV
!      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                            IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                           N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                           LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                              N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
                              A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
                              LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
END SUBROUTINE DGETRF

SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER*1        :: TRANS
      INTEGER            :: INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(FPK)      ::   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL(FPK)      :: array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) REAL(FPK)      :: array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)      ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL           ::  NOTRAN
!     ..
!     .. External Functions ..
!      LOGICAL            :: LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                     NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                     ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                     A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
END SUBROUTINE DGETRS

!

SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(FPK)      ::   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL(FPK)      :: array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      !    !    !    +    +    +       !    !    !   u14  u25  u36
!      !    !    +    +    +    +       !    !   u13  u24  u35  u46
!      !   a12  a23  a34  a45  a56      !   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   !      m21  m32  m43  m54  m65   !
!     a31  a42  a53  a64   !    !      m31  m42  m53  m64   !    !
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)      ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER          ::    I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
!      INTEGER          ::    IDAMAX
!      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 ) &
               CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
                           AB( KV+1, J ), LDAB-1 )
!
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J ) &
                  CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
                             AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
                             LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 )  INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of DGBTF2
!
END SUBROUTINE DGBTF2


SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER         :: INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(FPK)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL(FPK)      :: array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER           :: I, IP, IX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.EQ.0 ) RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I ) &
               CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I ) &
               CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
END SUBROUTINE DLASWP


INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    :: NAME, OPTS
      INTEGER            :: ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPE! for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPE!   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL

!  COmpiler dummies (R. Spurr)

      integer            nn3
      character*(1)      aa1
!     ..
!     .. Executable Statements ..
!

      nn3 = n3
      aa1 = opts(1:1)

!     Formerly a listed Goto: 
!       GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
      IF ( ISPEC.EQ.1 .or. ISPEC.eq.2 .or. ISPEC.eq.3 ) THEN

!  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) )    &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )  &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )  RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
!     Formerly a listed Goto:  GO TO ( 110, 200, 300 ) ISPEC
!
!  110 CONTINUE
      IF ( ISPEC .EQ. 1 ) THEN
!
!     ISPE! = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
!  200 CONTINUE
      ELSE IF ( ISPEC .EQ. 2 ) THEN
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.  &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
!  300 CONTINUE
      ELSE IF ( ISPEC .EQ. 3 ) THEN
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.  &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.  &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.  &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
!     END clause ISPEC = 1, 2 or 3
!
      ENDIF
!
!  400 CONTINUE
      ELSE IF ( ISPEC .EQ. 4 ) THEN
!
!     ISPE! = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
!  500 CONTINUE
      ELSE IF ( ISPEC .EQ. 5 ) THEN
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
!  600 CONTINUE 
      ELSE IF ( ISPEC .EQ. 6 ) THEN
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
!  700 CONTINUE
      ELSE IF ( ISPEC .EQ. 7 ) THEN
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
!  800 CONTINUE
      ELSE IF ( ISPEC .EQ. 8 ) THEN
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN

!
      ELSE
!
!     Invalid value of ISPEC
!
      ILAENV = -1
      RETURN
!
!     End of ILAENV
!
      ENDIF

!
!     End of ILAENV
!
END FUNCTION ILAENV

!

LOGICAL FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER         ::   INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME )  RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDI! machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDI! is assumed - ZCODE is the EBCDI! code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.  &
             INTA.GE.145 .AND. INTA.LE.153 .OR.  &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.  &
             INTB.GE.145 .AND. INTB.LE.153 .OR.  &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
END FUNCTION LSAME

SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        :: SRNAME
      INTEGER            :: INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
END SUBROUTINE XERBLA

!

SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(FPK)      ::   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(FPK)      :: array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(FPK)      ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER         ::   J, JP
!     ..
!     .. External Functions ..
!      INTEGER         ::   IDAMAX
!      EXTERNAL             IDAMAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )  &
               CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )  &
               CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,  &
                       A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
END SUBROUTINE DGETF2

SUBROUTINE dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(FPK)     ::  dx(*),dy(*)
      integer       :: i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return

      if(incx.eq.1.and.incy.eq.1)then

        m = mod(n,7)
        if( m .ne. 0 ) then
          do 30 i = 1,m
            dy(i) = dx(i)
   30     continue
          if( n .lt. 7 ) return
        endif

        mp1 = m + 1
        do 50 i = mp1,n,7
          dy(i) = dx(i)
          dy(i + 1) = dx(i + 1)
          dy(i + 2) = dx(i + 2)
          dy(i + 3) = dx(i + 3)
          dy(i + 4) = dx(i + 4)
          dy(i + 5) = dx(i + 5)
          dy(i + 6) = dx(i + 6)
   50   continue

      else
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
   10   continue

      endif

!  Finish

      return
END SUBROUTINE dcopy

SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1    ::   TRANSA, TRANSB
      INTEGER        ::   M, N, K, LDA, LDB, LDC
      REAL(FPK)      ::   ALPHA, BETA
!     .. Array Arguments ..
      REAL(FPK)      ::   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     ! := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and ! are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  ! an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(FPK)      :: array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL(FPK)      :: array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then ! need not be set on input.
!           Unchanged on exit.
!
!  !      - REAL(FPK)      :: array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  ! must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case ! need not be set on entry.
!           On exit, the array  !  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*! ).
!
!  LD!    - INTEGER.
!           On entry, LD! specifies the first dimension of ! as declared
!           in  the  calling  (sub)  program.   LD!  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL         :: LSAME
!      EXTERNAL           LSAME
!     .. External Subroutines ..
!      EXTERNAL           XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL        ::   NOTA, NOTB
      INTEGER        ::   I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(FPK)      ::   TEMP
!     .. Parameters ..
      REAL(FPK)      ::   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.   &
               ( .NOT.LSAME( TRANSA, 'C' ) ).AND.   &
               ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.   &
               ( .NOT.LSAME( TRANSB, 'C' ) ).AND.   &
               ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.   &
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )   &
         RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  ! := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  ! := alpha!A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  ! := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  ! := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
END SUBROUTINE DGEMM

SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,   &
                         BETA, Y, INCY )
!     .. Scalar Arguments ..
      REAL(FPK)      ::   ALPHA, BETA
      INTEGER        ::   INCX, INCY, LDA, M, N
      CHARACTER*1    ::   TRANS
!     .. Array Arguments ..
      REAL(FPK)      ::   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(FPK)      :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL(FPK)      :: array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL(FPK)      :: array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(FPK)      ::   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      REAL(FPK)      ::   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY

!     .. External Functions ..
!      LOGICAL        ::  LSAME
!      EXTERNAL           LSAME
!     .. External Subroutines ..
!      EXTERNAL           XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.   &
               .NOT.LSAME( TRANS, 'T' ).AND.   &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.   &
          ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )   &
         RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
END SUBROUTINE DGEMV

!

SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      REAL(FPK)      ::   ALPHA
      INTEGER        ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      REAL(FPK)      ::   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(FPK)      :: array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(FPK)      :: array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(FPK)      :: array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(FPK)   , PARAMETER   :: ZERO = 0.0D+0
!     .. Local Scalars ..
      REAL(FPK)      ::  TEMP
      INTEGER        ::  I, INFO, IX, J, JY, KX

!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )  RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
END SUBROUTINE DGER

!

SUBROUTINE dscal(n,da,dx,incx)

!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return

!        code for increment not equal to 1

      if(incx.ne.1) then
        nincx = n*incx
        do 10 i = 1,nincx,incx
          dx(i) = da*dx(i)
   10   continue
        return
      endif

!  clean-up loop

      m = mod(n,5)
      if( m .ne. 0 ) then
        do 30 i = 1,m
          dx(i) = da*dx(i)
   30   continue
        if( n .lt. 5 ) return
      endif

      mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return

END SUBROUTINE dscal

SUBROUTINE dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return

!       code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

        m = mod(n,3)
        if ( m .ne. 0 ) then
          do 30 i = 1,m
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
   30     continue
          if( n .lt. 3 ) return
        endif

        mp1 = m + 1
        do 50 i = mp1,n,3
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
          dtemp = dx(i + 1)
          dx(i + 1) = dy(i + 1)
          dy(i + 1) = dtemp
          dtemp = dx(i + 2)
          dx(i + 2) = dy(i + 2)
          dy(i + 2) = dtemp
   50   continue

      else
!
!       code for unequal increments or equal increments not equal to 1
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
          dtemp = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix = ix + incx
          iy = iy + incy
   10   continue
        return

      endif

      return
END SUBROUTINE dswap

!

SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
!     .. Scalar Arguments ..
      INTEGER            :: INCX, K, LDA, N
      CHARACTER*1        :: DIAG, TRANS, UPLO
!     .. Array Arguments ..
      REAL(FPK)          ::   A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - REAL(FPK)      :: array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - REAL(FPK)      :: array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(FPK)   , parameter   ::  ZERO = 0.0D+0
!     .. Local Scalars ..
      REAL(FPK)      ::   TEMP
      INTEGER        ::       I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL        ::       NOUNIT

!     .. External Functions ..
!      LOGICAL         ::      LSAME
!      EXTERNAL           LSAME
!     .. External Subroutines ..
!      EXTERNAL           XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
               .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND. &
               .NOT.LSAME( TRANS, 'T' ).AND. &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND. &
               .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTBSV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )  RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT ) X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT ) X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := inv( A')*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )  TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )  TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )  KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT ) TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )  KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DTBSV  .
!
END SUBROUTINE DTBSV

!

SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )

!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL(FPK)      ::   ALPHA
!     .. Array Arguments ..
      REAL(FPK)      ::   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL(FPK)      :: array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL(FPK)      :: array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL            :: LSAME
!      EXTERNAL              LSAME
!     .. External Subroutines ..
!      EXTERNAL              XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC             MAX
!     .. Local Scalars ..
      LOGICAL            :: LSIDE, NOUNIT, UPPER
      INTEGER            :: I, INFO, J, K, NROWA
      REAL(FPK)          ::   TEMP
!     .. Parameters ..
      REAL(FPK)   , PARAMETER  ::  ONE  = 1.0D+0
      REAL(FPK)   , PARAMETER  ::  ZERO= 0.0D+0
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND. &
               ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND. &
               ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
               ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )  B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT ) TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!    End of DTRSM .
! 
END SUBROUTINE DTRSM

INTEGER FUNCTION idamax(n,dx,incx)
!
      implicit none

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision :: dx(*),dmax
      integer          :: i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return

!   code for increment equal to 1

      if(incx.eq.1)then

        dmax = dabs(dx(1))
        do 30 i = 2,n
          if(dabs(dx(i)) .gt. dmax) then
            idamax = i
            dmax = dabs(dx(i))
          endif
   30   continue

!  code for increment not equal to 1

      else
        ix = 1
        dmax = dabs(dx(1))
        ix = ix + incx
        do 10 i = 2,n
           if(dabs(dx(ix)).gt.dmax) then
             idamax = i
             dmax = dabs(dx(ix))
           endif
          ix = ix + incx
   10   continue

      endif

      return
END FUNCTION idamax

!

LOGICAL FUNCTION XFINDPAR ( NUNIT, PARNAME )

!  FINDPAR input subroutines.

!----- The basic idea

!  FINDPAR is a poor man's RDPAR.  It implements some of the
! central features of RDPAR --- the ability to find a named set
! of data values in an input file, the relaxation of the need for
! data in the input file to appear in the order in which it's
! used, the provision for comments in the input file --- in
! (almost) standard FORTRAN 77.

!  Among the features of RDPAR that are not implemented in FINDPAR:
! all redirection of input (from one parameter name to another, from
! file to terminal, from file to command line); freer-format data
! specification on input; tracing options.

!----- Using FINDPAR

!  To use FINDPAR, you open your input file using ordinary FORTRAN
! procedures.  You will also use ordinary FORTRAN to read data
! values.  The only difference is that you can use the FINDPAR call
! at any point to position the file to a point after a line containing
! a parameter name.  These parameter names help to document what the
! input values are for; they also don't have to appear in the order in
! which the program uses them.  An input file might look like this:

! ! Parameters for July 16 run.

! NLAYERS
! 5

! LAYER PARAMETERS
! 5,1013,296
! .001,.001,.001,.001,.001,.001,.001,.001

! MOLECULE FLAGS
! t,t,t,t,t,t,t,t

!  By convention, ! or a space introduces a comment line.  A
! parameter name appears on a line of its own, followed by the data
! lines to be read using FORTRAN statements.  Here is a sample of
! FORTRAN code to read this parameter file:

!      OPEN (UNIT = IN_UNIT, FILE = INPUT, STATUS = 'OLD')
!      CALL FINDPAR ( IN_UNIT, 'NLAYERS' )
!      READ (IN_UNIT, *) NLAYERS
!      CALL FINDPAR ( IN_UNIT, 'LAYER PARAMETERS' )
!      READ (IN_UNIT, *) Z, PRESS, TEMP
!      READ (IN_UNIT, *) BRO, CLO, HCHO, NO2, O2, O3, OCLO, SO2
!      CALL FINDPAR ( IN_UNIT, 'MOLECULE FLAGS' )
!      READ (IN_UNIT, *) IF_BRO, IF_CLO, IF_HCHO, IF_NO2, IF_O2, IF_O3,
!      IF_OCLO, IF_SO2
!      CLOSE (UNIT = IN_UNIT)

!----- Details of parameter file format

! The parameter file consists of a sequence of comments and
! parameter-name/data chunks, in any order (except as noted below).

! Comments are essentially anything that doesn't happen to get matched
! as a parameter name or read as data.  However, by convention comments
! are introduced by an exclamation point, to make absolutely sure that
! they don't get mistaken for parameter names.  It is also conventional
! to ``comment out'' a parameter name that you don't want to match by
! simply indenting it one space.

! Parameter names can theoretically be almost any sequence of
! characters.  However, FINDPAR ignores trailing spaces both in the
! parameter name you give it and in lines it reads from the parameter
! file.  And, by convention, parameter names don't begin with spaces.
! Parameter names may contain embedded spaces.  Parameter name matching
! is case-insensitive but otherwise exact: if you use odd things
! like control characters in a parameter name they have to match
! exactly, too, so it's usually best to avoid control characters and
! tabs in parameter names.

! The data lines that follow a line with a parameter name are ordinary
! FORTRAN input lines.  FINDPAR has no influence on how they're handled,
! and all the ordinary rules about FORTRAN READ statements apply.

! Comments may not appear on the same lines as parameter names or data,
! and they can't appear within parameter name/data chunks.  When FINDPAR
! locates a parameter name it positions the file at the line following
! that name; your program will then use ordinary FORTRAN READ statements
! to read the data, and will probably be confused by comments that appear
! at that point.

! There's a maximum length for input lines and parameter names, set
! by the parameter MAXLINE in the code.

! By convention, a parameter name that ends in ? precedes an
! input line that contains a single logical variable (T or F).

!----- Repeated parameters

! Sometimes you want to read several parameters with the same name:
! for example, parameters for each layer of the atmosphere, where there
! might be a large number of layers.  In this case the order of
! parameters in the file does matter, since you usually just want to
! put the numbers in the file in some natural order rather than giving
! a number in each input line specifying where it appears in the input
! order.

! The normal FINDPAR approach is to organize these parameters as
! follows: begin with some parameter that appears only once in the file
! (the number of layers, for example); then follow it with the instances
! of the repeated parameter and associated data.  The only-once
! parameter can even be a dummy, with no data following; its importance
! is that reading it gets the parameter file positioned at the right
! point.

! The other problem here is knowing when a list of repeated parameters
! is over.  FINDPAR doesn't have any of the tricks for doing this that
! RDPAR has; you must either know exactly how many instances of the
! repeated parameter there are going to be (by reading some once-only
! parameter that tells you), or else you need some special numbers to
! flag the end of the list (zeros, negative numbers, etc.).  You can't
! just keep looking for the same parameter name indefinitely, because
! FINDPAR will just rewind the file for you and loop through its
! contents forever.

!----- Errors and XFINDPAR

! If FINDPAR can't find a parameter, it prints an error message
! and returns.  Your READ statements following the FINDPAR call are
! very likely to run into errors of their own at this point, since
! you'll be positioned at some random point in the file (usually at
! the beginning, in this version, but that isn't true in all cases).

! If you want to do something more intelligent about such errors,
! you can use XFINDPAR, which is a function returning a logical value:
! .true. if the parameter name was found, .false. if not; XFINDPAR
! doesn't display FINDPAR's error message when a parameter name cannot
! be found.

! Both FINDPAR and XFINDPAR can run into ordinary FORTRAN input errors
! as they read the file: no special action is taken on these---the
! system's default action, whatever that is, occurs.

!       9/10/91         John Lavagnino

!  Find a parameter name in a parameter file.

!        SUBROUTINE FINDPAR ( NUNIT, PARNAME )

!  Input arguments.

!        INTEGER         NUNIT
!        CHARACTER * (*) PARNAME

!  Local variables.

!        LOGICAL         XFINDPAR
!        EXTERNAL        XFINDPAR

!  No local variables need to be SAVEd.

!***********************************************************************

!        IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
!          WRITE (*, *) 'FINDPAR error'
!          WRITE (*, *) '   Unit number ', NUNIT
!          WRITE (*, *) '   Parameter name ', PARNAME
!        END IF

!        END
!
!  Find a parameter name in a paramter file: return .true. if found,
! .false. if not.

        implicit none

!  Parameters.

        INTEGER, parameter  :: MAXLINE = 132

!  Input arguments.

        INTEGER           :: NUNIT
        CHARACTER (LEN=*) :: PARNAME

!  Local variables.

        INTEGER                 :: III
        INTEGER                 :: PARLEN
        INTEGER                 :: START
        INTEGER                 :: END
        LOGICAL                 :: ENDSEEN
        CHARACTER * (MAXLINE)   :: LINE
        CHARACTER * (MAXLINE)   :: NAME

!  No local variables need to be SAVEd.

!***********************************************************************

!  Determine the length of the parameter name.

        DO 10 III = LEN ( PARNAME ), 1, -1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            GO TO 20
          END IF
10      CONTINUE

!  If we get to here, then name contains nothing but blanks.
!  We just return, claiming success, in such a case.

        GO TO 500

!  If we get here, then there's a non-null parameter name; but it
!  might still be too long, in which case we always return
!  signaling failure, since we couldn't ever find such a name in the
!  file.

20      CONTINUE
        IF ( PARLEN  .GT.  MAXLINE ) GO TO 400

!  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

!  Top of main loop.

        ENDSEEN = .FALSE.

100     CONTINUE

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )
          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) ) GO TO 100
          START = PARLEN + 1
          END = MAXLINE
          IF ( START  .GT.  END ) GO TO 500
          IF ( LINE ( START : END ) .EQ. ' ' ) GO TO 500
          GO TO 100

!  End-of-file branch.

200       CONTINUE
          REWIND ( UNIT = NUNIT )
          IF ( ENDSEEN ) THEN
            GO TO 400
          ELSE
            ENDSEEN = .TRUE.
            GO TO 100
          END IF

!  End of loop: failure, no parameter name found.

400     CONTINUE
        XFINDPAR = .FALSE.
        RETURN

!  End of loop: successful location of parameter name.

500     CONTINUE
        XFINDPAR = .TRUE.
        RETURN

END FUNCTION XFINDPAR

!

LOGICAL FUNCTION XFINDPAR_NOG ( NUNIT, PARNAME )

!  Find a parameter name in a parameter file: return .true. if found,
! .false. if not.

        implicit none

!  Parameters.

        INTEGER, parameter  :: MAXLINE = 132

!  Input arguments.

        INTEGER         :: NUNIT
        CHARACTER * (*) :: PARNAME

!  Local variables.

        INTEGER                 :: III
        INTEGER                 :: PARLEN
        INTEGER                 :: START
        INTEGER                 :: ENDL
        LOGICAL                 :: ENDSEEN
        LOGICAL                 :: LOOP
        LOGICAL                 :: CARRYON
        CHARACTER * (MAXLINE)   :: LINE
        CHARACTER * (MAXLINE)   :: NAME

!  No local variables need to be SAVEd.

!***********************************************************************

!  Determine the length of the parameter name.

        LOOP = .TRUE.
        III = LEN ( PARNAME ) + 1
        DO WHILE (LOOP.AND.III.GT.1)
          III = III - 1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            LOOP = .false.
          END IF
        ENDDO

!  If we get to here, then name contains nothing but blanks.
!  We just return, claiming success, in such a case.

        IF ( LOOP ) THEN
          XFINDPAR_NOG = .true.
          RETURN
        ENDIF

!  If we get here, then there's a non-null parameter name; but it
!  might still be too long, in which case we always return
!  signaling failure, since we couldn't ever find such a name in the
!  file. Bo parameter name found.

        IF ( PARLEN  .GT.  MAXLINE ) THEN
          XFINDPAR_NOG = .FALSE.
          RETURN
        ENDIF

!  If we get here, then there's a non-null parameter name; and
!  It is not too long, so now examine it.

!  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

!  Top of main loop.

        ENDSEEN = .FALSE.
        CARRYON = .false.

        DO WHILE ( .not. ENDSEEN )

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )

          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) ) THEN
            START = 0
            CARRYON = .true.
          ELSE
            START = PARLEN + 1
            ENDL = MAXLINE
            IF ( START  .GT.  ENDL )  THEN
               XFINDPAR_NOG = .true.
               RETURN
            ENDIF
            IF ( LINE ( START : ENDL ) .EQ. ' ' ) THEN
               XFINDPAR_NOG = .true.
               RETURN
            ENDIF
          ENDIF

          CARRYON = .true.

!  End-of-file branch.

 200      continue
          IF ( .not. CARRYON ) THEN
            REWIND ( UNIT = NUNIT )
            IF ( ENDSEEN ) THEN
               XFINDPAR_NOG = .FALSE.
               RETURN
            ELSE
              ENDSEEN = .TRUE.
            ENDIF
          END IF

!  End of loop

        ENDDO

        RETURN 
END FUNCTION XFINDPAR_NOG

!

LOGICAL FUNCTION GFINDPAR (NUNIT, PREFIX, ERROR, PARNAME)

!  Find a parameter name in a parameter file, with optional added
! prefix on parameter name.  Error messages are displayed but
! instead of stopping the program a logical flag is returned.
!  Returns .TRUE. if some form of the parameter name was found;
! .FALSE. if not.
! ERROR and GFINDPAR are initialised for a successful search.
! ERROR is just the negation of GFINDPAR. Error messages are removed. 

        implicit none

!  Input arguments.  If PREFIX is something other than a bunch of
! blanks, it will be added onto PARNAME, with a blank to separate
! it, for the FINDPAR call.  If that FINDPAR call fails, then
! PARNAME without prefix is tried instead, and only if that call fails
! is there an error message.  This means that unprefixed versions of
! the parameter name act as defaults.
!  If PREFIX is a bunch of blanks, this works much as FINDPAR does.
!  (This is assuming that PREFIX has no trailing blanks.  Probably
! it should be checking for them and trimming them if necessary.)

        INTEGER           :: NUNIT
        CHARACTER (LEN=*) :: PREFIX
        CHARACTER (LEN=*) :: PARNAME

!  Modified argument.  If all goes well, this is not changed.  If
! the parameter name can't be found, it's set to .TRUE.

        LOGICAL         :: ERROR

!  Local variables.

        CHARACTER * 132 :: LINE

!***********************************************************************

!  Initialise return values.

        GFINDPAR = .TRUE.
        ERROR    = .FALSE.

!  Compose search line for first search.

        IF ( PREFIX  .EQ.  ' ' ) THEN
          LINE = PARNAME
        ELSE
          LINE = PREFIX // ' ' // PARNAME
        END IF

!  Try the parameter name with prefix.

        IF (XFINDPAR (NUNIT, LINE)) THEN
          RETURN
        END IF

!  If that doesn't work, try it without the prefix, if the prefix
! was non-null.

        IF (PREFIX  .NE.  ' ') THEN
          IF (XFINDPAR (NUNIT, PARNAME)) THEN
            RETURN
          END IF
        END IF

!  If we got here, we just couldn't find the parameter name.

        GFINDPAR = .FALSE.
        ERROR    = .TRUE.

END FUNCTION GFINDPAR

!

SUBROUTINE FINDPAR_ERROR (  ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_MESSAGES, LIDORT_SERIOUS

      implicit none

!  subroutine input arguments

      LOGICAL, intent(in)       :: ERROR
      CHARACTER (LEN=*), intent(in) :: PAR_STR
      
!  subroutine Output or In/Out arguments

      INTEGER, intent(inout)       :: STATUS
      INTEGER, intent(inout)       :: NM
      CHARACTER (LEN=*), intent(inout) :: MESSAGES(0:MAX_MESSAGES), ACTIONS(0:MAX_MESSAGES)

      IF ( ERROR ) THEN
        NM = NM + 1
        STATUS = LIDORT_SERIOUS
        MESSAGES(NM) = 'Cannot find string: '//PAR_STR(1:LEN_STRING(PAR_STR))
        ACTIONS(NM)  = 'Check Spelling of String in input file'
      ENDIF

!  finish

      RETURN
END SUBROUTINE FINDPAR_ERROR

!  End

end module lidort_aux_m
