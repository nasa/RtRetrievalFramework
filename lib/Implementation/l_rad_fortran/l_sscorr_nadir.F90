module l_sscorr_nadir_m

!  For given wavelength, routine calculates TOA First Order upwelling
!  Stokes vectors and any number of LPLS Jacobians (profile/surface)
!  for the Atmospheric Single-scatter and Surface Direct-Beam solar sources

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Regular-PS calculations

!  This is based on Rob's standalone FO code
!  Modified to remove thermal, downwelling and plane-parallel/enhanced ps
!  SS and DB output combined

!  17 April 2013, V. Natraj, JPL

!  The subroutines are
!       l_sscorr_nadir

!  All subroutines public

public

contains

subroutine l_sscorr_nadir &
   ( nlayers, nstokes,                              & ! Inputs (Layers/Stokes control)
     linearize, s_linearize,                        & ! Inputs (Linearization flags)
     n_wfuncs, n_s_wfuncs,                          & ! Inputs (dimensioning)
     do_lambertian,                                 & ! Inputs (Flags - General)
     r1_surf, zmat_up, flux,                        & ! Inputs (Optical)
     deltau_vert_i, omega_total_i,                  & ! Inputs (Optical - Intensity)
     deltau_vert, omega_total,                      & ! Inputs (Optical - Polarization)
     Ls_r1_surf, L_zmat_up,                         & ! Inputs (Optical - Linearized)
     L_deltau_vert_i, L_omega_total_i,              & ! Inputs (Optical - Linearized - Intensity)
     L_deltau_vert, L_omega_total,                  & ! Inputs (Optical - Linearized - Polarization)
     Mu0, Mu1, height_grid,                         & ! Inputs (Geometry)
     sunpaths, ntraverse,                           & ! Inputs (Geometry)
     Stokes_SS, L_Stokes_SS, LS_Stokes_SS )                              ! Outputs

!  Stand alone routine for TOA SS field (Radiances and Jacobians) with Solar sources alone

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: n_wfuncs, n_s_wfuncs

!  General flags

      LOGICAL, Intent(in) ::  DO_LAMBERTIAN

!  Jacobian Flags

      LOGICAL, Intent(in) ::  linearize
      LOGICAL, Intent(in) ::  s_linearize

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: DELTAU_VERT_I  ( NLAYERS )
      REAL(fpk), Intent(in) :: OMEGA_TOTAL_I  ( NLAYERS )
      REAL(fpk), Intent(in) :: DELTAU_VERT    ( NLAYERS )
      REAL(fpk), Intent(in) :: OMEGA_TOTAL    ( NLAYERS )
      REAL(fpk), Intent(in) :: ZMAT_UP        ( NLAYERS, NSTOKES )

!  Surface reflectivity (Could be the albedo)

      REAL(fpk), Intent(in) :: R1_SURF(nstokes), FLUX

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_DELTAU_VERT_I    ( NLAYERS, n_wfuncs )
      REAL(fpk), Intent(in) :: L_OMEGA_TOTAL_I    ( NLAYERS, n_wfuncs )
      REAL(fpk), Intent(in) :: L_DELTAU_VERT      ( NLAYERS, n_wfuncs )
      REAL(fpk), Intent(in) :: L_OMEGA_TOTAL      ( NLAYERS, n_wfuncs )
      REAL(fpk), Intent(in) :: L_ZMAT_UP          ( NLAYERS, NSTOKES, n_wfuncs )
      REAL(fpk), Intent(in) :: LS_R1_SURF         ( NSTOKES, n_s_wfuncs )

!  Geometrical inputs
!  ------------------

!  Mu0 = cos(theta_boa), required for surface term
!  Mu1 = cos(alpha_boa)
!  Height grid

      real(fpk), Intent(in)  :: Mu0
      real(fpk), Intent(in)  :: Mu1
      real(fpk), Intent(in)  :: height_grid(0:nlayers)

!  solar paths

      integer,   Intent(in)  :: ntraverse  (0:nlayers)
      real(fpk), Intent(in)  :: sunpaths   (0:nlayers,nlayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: Stokes_SS    ( nstokes )
      real(fpk), Intent(Out)  :: L_Stokes_SS  ( nstokes, nlayers, n_wfuncs )
      real(fpk), Intent(Out)  :: LS_Stokes_SS ( nstokes, n_s_wfuncs )

!  LOCAL
!  -----

!  Extinctions

      REAL(fpk)  :: EXTINCTION_I   ( NLAYERS )
      REAL(fpk)  :: EXTINCTION     ( NLAYERS )
      REAL(fpk)  :: L_EXTINCTION   ( NLAYERS, n_wfuncs )
      REAL(fpk)  :: L_EXTINCTION_I ( NLAYERS, n_wfuncs )

!  Attenuations

      real(fpk)  :: suntau_i              (0:nlayers)
      real(fpk)  :: attenuations_i        (0:nlayers)

      real(fpk)  :: suntau                (0:nlayers)
      real(fpk)  :: attenuations          (0:nlayers)

      real(fpk)  :: L_suntau_i            (0:nlayers,0:nlayers,n_wfuncs)
      real(fpk)  :: L_attenuations_i      (0:nlayers,0:nlayers,n_wfuncs)

      real(fpk)  :: L_suntau              (0:nlayers,0:nlayers,n_wfuncs)
      real(fpk)  :: L_attenuations        (0:nlayers,0:nlayers,n_wfuncs)

!  Source function integration results

      real(fpk)  :: sources_up       ( nlayers, nstokes )
      real(fpk)  :: lostrans_up_i    ( nlayers )
      real(fpk)  :: multiplier_up_i  ( nlayers )
      real(fpk)  :: lostrans_up      ( nlayers )
      real(fpk)  :: multiplier_up    ( nlayers )

      real(fpk)  :: L_sources_up      ( nlayers,nstokes,0:nlayers,n_wfuncs )
      real(fpk)  :: L_lostrans_up_i   ( nlayers, n_wfuncs )
      real(fpk)  :: L_multiplier_up_i ( nlayers,0:nlayers,n_wfuncs )
      real(fpk)  :: L_lostrans_up     ( nlayers, n_wfuncs )
      real(fpk)  :: L_multiplier_up   ( nlayers,0:nlayers,n_wfuncs )

!  Local cumulative source terms

      real(fpk)  :: cumsource_db      ( 0:nlayers, nstokes )
      real(fpk)  :: cumsource_up      ( 0:nlayers, nstokes )

      real(fpk)  :: L_cumsource       ( nstokes, n_wfuncs )
      real(fpk)  :: LS_cumsource      ( nstokes, n_s_wfuncs )

!  Help

      integer    :: n, ns, k, q, o1, nstart, nc

      real(fpk)  :: help, factor1, factor2, factor3, m4, m4a, rhelp(4), shelp(4)
      real(fpk)  :: factor1_i, factor2_i, factor3_i
      real(fpk)  :: L_factor1, L_factor2
      real(fpk)  :: L_factor1_i, L_factor2_i
      real(fpk)  :: sumd_i
      real(fpk)  :: sumd
      real(fpk)  :: lostau, L_lostau(n_wfuncs)
      real(fpk)  :: lostau_i, L_lostau_i(n_wfuncs)
      real(fpk)  :: L_Shelp

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      STOKES_SS       = zero
      L_STOKES_SS     = zero
      LS_STOKES_SS    = zero

!  Zero local sources

      sources_up   = zero
      L_sources_up = zero
      lostrans_up_i = zero   ; multiplier_up_i   = zero
      lostrans_up   = zero   ; multiplier_up   = zero
      L_lostrans_up_i = zero ; L_multiplier_up_i = zero
      L_lostrans_up = zero   ; L_multiplier_up = zero

!  Bookkeeping

      ns = nstokes

!  Use basic definitions

      do n = 1, nlayers
         help = height_grid(n-1) - height_grid(n)
         extinction(n) = deltau_vert(nlayers-n+1) / help
         extinction_i(n) = deltau_vert_i(nlayers-n+1) / help
      enddo

!  Linearized extinctions

      if ( linearize ) THEN   
         do k = 1, nlayers
            help = height_grid(k-1) - height_grid(k)
            do q = 1, n_wfuncs
               l_extinction(k,q) = l_deltau_vert(nlayers-k+1,q) / help
               l_extinction_i(k,q) = l_deltau_vert_i(nlayers-k+1,q) / help
            enddo
         enddo
      endif

!  Attenuations
!  ============

!  Initialize

      Attenuations     = zero
      Attenuations_i   = zero
      L_Attenuations   = zero
      L_Attenuations_i = zero
      Suntau           = zero ; L_suntau              = zero
      Suntau_i         = zero ; L_suntau_i            = zero
      nstart = nlayers

!  Attenuations to End points (including TOA). All representations
!  MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         sumd_i = ZERO
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
            sumd_i = sumd_i + extinction_i(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         suntau_i(n) = sumd_i
         if (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if (sumd_i .lt. cutoff ) Attenuations_i(n) = exp( - sumd_i )
         if ( linearize ) then
            do k = 1, nlayers
               if ( k.le.ntraverse(n) ) then
                  do q = 1, n_wfuncs
                     L_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k)
                     L_suntau_i(n,k,q) = L_extinction_i(k,q) * sunpaths(n,k)
                     L_Attenuations(n,k,q) = - Attenuations(n) * L_suntau(n,k,q)
                     L_Attenuations_i(n,k,q) = - Attenuations_i(n) * L_suntau_i(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  Adjust nstart

      do n = 1, nlayers
         if ( attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Layer integrated Solar sources
!  ==============================

!  Regular PS multipliers and LOSTRANS
!  -----------------------------------

      do n = nlayers, 1, -1

!  LOS transmittance (not for the Horizontal case)

         if ( Mu1 .gt. zero ) then
            lostau = deltau_vert(nlayers-n+1)/Mu1
            lostau_i = deltau_vert_i(nlayers-n+1)/Mu1
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( lostau_i .lt. cutoff ) lostrans_up_i(n) = exp( - lostau_i )
            if ( linearize ) then
               do q = 1, n_wfuncs
                  L_lostau(q) = L_deltau_vert(nlayers-n+1,q) / Mu1
                  L_lostau_i(q) = L_deltau_vert_i(nlayers-n+1,q) / Mu1
                  L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  L_lostrans_up_i(n,q) = - L_lostau_i(q) * lostrans_up_i(n)
               enddo
            endif
         endif

!  Multipliers

         if ( n.le.nstart  ) then
            if ( Mu1 .gt. zero ) then
               if ( attenuations(n-1).ne.zero ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  factor1_i = Attenuations_i(n-1) - Attenuations_i(n)*lostrans_up_i(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  factor2_i = (suntau_i(n) - suntau_i(n-1))/lostau_i
                  factor3 = one + factor2
                  factor3_i = one + factor2_i
                  multiplier_up(n) = factor1 / factor3
                  multiplier_up_i(n) = factor1_i / factor3_i
                  if ( linearize ) then
                     do k = 1, nlayers
                        if ( k.le.ntraverse(n) ) then
                           do q = 1, n_wfuncs
                              L_factor1 = L_Attenuations(n-1,k,q) - L_Attenuations(n,k,q)*lostrans_up(n)
                              L_factor1_i = L_Attenuations_i(n-1,k,q) - L_Attenuations_i(n,k,q)*lostrans_up_i(n)
                              L_factor2 = ( L_suntau(n,k,q) - L_suntau(n-1,k,q) ) / lostau
                              L_factor2_i = ( L_suntau_i(n,k,q) - L_suntau_i(n-1,k,q) ) / lostau_i
                              if ( k.eq.n ) then
                                 L_factor1 = L_factor1 - Attenuations(n) * L_lostrans_up(n,q)
                                 L_factor1_i = L_factor1_i - Attenuations_i(n) * L_lostrans_up_i(n,q)
                                 L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                                 L_factor2_i = L_factor2_i - factor2_i *  L_lostau_i(q) / lostau_i
                              endif
                              L_multiplier_up(n,k,q) = ( L_factor1 - multiplier_up(n) * L_factor2 ) / factor3
                              L_multiplier_up_i(n,k,q) = ( L_factor1_i - multiplier_up_i(n) * L_factor2_i ) / factor3_i
                           enddo
                        endif
                     enddo
                  endif
               endif
            else
               if ( attenuations(n-1).ne.zero ) then
                  multiplier_up(n) = Attenuations(n-1)
                  multiplier_up_i(n) = Attenuations_i(n-1)
                  if ( linearize ) then
                     do k = 1, nlayers
                        if ( k.le.ntraverse(n) ) then
                           do q = 1, n_wfuncs
                              L_multiplier_up(n,k,q) = L_Attenuations(n-1,k,q)
                              L_multiplier_up_i(n,k,q) = L_Attenuations_i(n-1,k,q)
                           enddo
                        endif
                     enddo
                  endif
               endif
            endif
         endif

!  End layer loop

      enddo

!  Layer sources
!  -------------

      do n = nlayers, 1, -1
         nc = nlayers - n + 1
         if ( n.le.nstart ) then
            do o1 = 1, 1
               shelp(o1) = zmat_up(nc,o1) * omega_total_i(nc)
               sources_up(n,o1) = shelp(o1) * multiplier_up_i(n)
            enddo
            do o1 = 2, nstokes
               shelp(o1) = zmat_up(nc,o1) * omega_total(nc)
               sources_up(n,o1) = shelp(o1) * multiplier_up(n)
            enddo
            if ( linearize ) then
               do k = 1, nlayers
                  do q = 1, n_wfuncs
                     do o1 = 1, 1
                        L_sources_up(n,o1,k,q) = shelp(o1) * L_multiplier_up_i(n,k,q)
                        if ( k.eq.n ) then
                           L_Shelp = L_zmat_up(nc,o1,q) * omega_total_i(nc) + zmat_up(nc,o1) * L_omega_total_i(nc,q)
                           L_sources_up(n,o1,k,q) =  L_sources_up(n,o1,k,q) + L_Shelp * multiplier_up_i(n)
                        endif
                     enddo
                     do o1 = 2, nstokes
                        L_sources_up(n,o1,k,q) = shelp(o1) * L_multiplier_up(n,k,q)
                        if ( k.eq.n ) then
                           L_Shelp = L_zmat_up(nc,o1,q) * omega_total(nc) + zmat_up(nc,o1) * L_omega_total(nc,q)
                           L_sources_up(n,o1,k,q) =  L_sources_up(n,o1,k,q) + L_Shelp * multiplier_up(n)
                        endif
                     enddo
                  enddo
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB(NC,:) = zero

!  Surface term

      RHELP = zero; M4 = 4.0_fpk * MU0 ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * R1_SURF(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations_i(nlayers)
      else
         do o1 = 1, 1
            RHELP(O1) = M4 * R1_SURF(O1)
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations_i(nlayers)
         enddo
         do o1 = 2, nstokes
            RHELP(O1) = M4 * R1_SURF(O1)
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
         enddo
      endif

!  Cumulative source terms : Loop over layers working upwards to TOA

      DO N = NLAYERS, 1, -1
         NC = NLAYERS + 1 - N
         do o1 = 1, 1
            CUMSOURCE_DB(NC,O1) = LOSTRANS_UP_I(N) * CUMSOURCE_DB(NC-1,O1)
            CUMSOURCE_UP(NC,O1) = LOSTRANS_UP_I(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
         enddo
         do o1 = 2, nstokes
            CUMSOURCE_DB(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1,O1)
            CUMSOURCE_UP(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
         enddo
      ENDDO
      do o1 = 1, nstokes
         STOKES_SS   (O1) = FLUX * (CUMSOURCE_UP(NC,O1) + CUMSOURCE_DB(NC,O1))
      enddo

!  Surface WFs

      if ( s_linearize ) then
         LS_CUMSOURCE = zero
         do q = 1, n_s_wfuncs
            if ( DO_LAMBERTIAN ) then
               LS_cumsource(1,q) = M4A * LS_R1_SURF(1,q)
            else
               do o1 = 1, nstokes
                  LS_cumsource(o1,q) = M4A * LS_R1_SURF(O1,q)
               enddo
            endif
         enddo
         DO N = NLAYERS, 1, -1
            do q = 1, n_s_wfuncs
               LS_cumsource(1,q) = LOSTRANS_UP_I(N) * LS_CUMSOURCE(1,q)
               LS_cumsource(2:ns,q) = LOSTRANS_UP(N) * LS_CUMSOURCE(2:ns,q)
            enddo
         ENDDO
         do q = 1, n_s_wfuncs
            LS_STOKES_SS(1:ns,Q) = FLUX * LS_CUMSOURCE(1:ns,Q)
         enddo
      endif

!  Profile Wfs (atmospheric term)

      if ( linearize ) then
         do k = 1, nlayers
            L_CUMSOURCE = zero
            DO N = NLAYERS, 1, -1
               NC = NLAYERS + 1 - N
               if ( k.eq.n ) then
                  do q = 1, n_wfuncs
                     L_cumsource(1,q) = L_SOURCES_UP(N,1,K,Q)    + &
                          L_LOSTRANS_UP_I(N,Q) * CUMSOURCE_UP(NC-1,1) + &
                            LOSTRANS_UP_I(N)   * L_CUMSOURCE(1,Q)
                     L_cumsource(2:ns,q) = L_SOURCES_UP(N,2:ns,K,Q)    + &
                          L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,2:ns) + &
                            LOSTRANS_UP(N)   * L_CUMSOURCE(2:ns,Q)
                  enddo
               else
                  do q = 1, n_wfuncs
                     L_cumsource(1,q) = L_SOURCES_UP(N,1,K,Q)  + &
                            LOSTRANS_UP_I(N)   * L_CUMSOURCE(1,Q)
                     L_cumsource(2:ns,q) = L_SOURCES_UP(N,2:ns,K,Q)  + &
                            LOSTRANS_UP(N)   * L_CUMSOURCE(2:ns,Q)
                  enddo
               endif
            ENDDO
            do q = 1, n_wfuncs
               L_STOKES_SS(1:ns,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
         enddo
      endif

!  Profile Wfs (direct beam term)

      if ( linearize ) then
         do k = 1, nlayers
            do q = 1, n_wfuncs
               L_CUMSOURCE(1,q) = RHELP(1) * L_attenuations_i(nlayers,k,q)
               L_CUMSOURCE(2:ns,q) = RHELP(2:ns) * L_attenuations(nlayers,k,q)
            enddo
            DO N = NLAYERS, 1, -1
               NC = NLAYERS + 1 - N
               if ( k.eq.n ) then
                  do q = 1, n_wfuncs
                     L_cumsource(1,q) =  L_LOSTRANS_UP_I(N,Q) * CUMSOURCE_DB(NC-1,1) + &
                                              LOSTRANS_UP_I(N)   * L_CUMSOURCE(1,Q)
                     L_cumsource(2:ns,q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1,2:ns) + &
                                              LOSTRANS_UP(N)   * L_CUMSOURCE(2:ns,Q)
                  enddo
               else
                  do q = 1, n_wfuncs
                     L_cumsource(1,q) = LOSTRANS_UP_I(N) * L_CUMSOURCE(1,Q)
                     L_cumsource(2:ns,q) = LOSTRANS_UP(N) * L_CUMSOURCE(2:ns,Q)
                  enddo
               endif
            ENDDO
            do q = 1, n_wfuncs
               L_STOKES_SS(1:ns,K,Q) = L_STOKES_SS(1:ns,K,Q) + FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
         enddo
      endif

!  Finish

      return
end subroutine l_sscorr_nadir

!  End module

end module l_sscorr_nadir_m
