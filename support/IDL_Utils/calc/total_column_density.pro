FUNCTION Total_Column_Density, pressure, temperature, profile, ALTITUDE=altitude, LATITUDE=latitude, PSURF=psurf, TOT_COL_AIR=tot_col_air

   avogadro = 6.02e23 ; kg * mole

   boltzmann =  1.3806503e10-23 ; m^2 kg s^-2 K^-1

   mass_p_molec_air = .02896 / avogadro  ; kg / molecule

   n_levels = n_elements(profile)

   ;; Interpolate mixing ratio to surface pressure
   if n_elements(psurf) ne 0 then begin
       profile[n_levels-1] =                                           $
         (profile[n_levels-1] - profile[n_levels-2])                   $
         / (alog(pressure[n_levels-1])                              $
            - alog(pressure[n_levels-2]))                            $
         * (alog(psurf) - alog(pressure[n_levels-2]))  $
         + profile[n_levels-2]
       pressure[n_levels-1] = psurf
   endif

   if n_elements(altitude) le 0 then begin
       Message, /CONTINUE, 'WARNING: Using 0.0 for altitude. This may be probably wrong for your calculation.'
       altitude = 0.0
   endif

   tot_col_species = 0.0E0
   tot_col_air = 0.0E0

   level_heights = geo_height(pressure[*], temperature[*], profile[*], altitude, latitude=latitude, gravity=gravity)

   if level_heights[0] EQ -1 then $
     Message, 'Geopotential function failed'

   totp = pressure[n_levels-1] - pressure[0]
   for lev_idx = 0, n_levels-2 do begin
       ;; Calculate interpolation factor for layer
       dpres = pressure[lev_idx + 1] - pressure[lev_idx]
       
       rp = pressure[lev_idx + 1] / pressure[lev_idx]
       interp_f = pressure[lev_idx + 1] / dpres - 1.0D / alog(rp)

       ;; Compute mean gravity for layer
       mean_gravity = gravity[lev_idx] + (gravity[lev_idx + 1] - gravity[lev_idx]) * interp_f

       ;; Calculate layer total number air density
       delta_pres = pressure[lev_idx]  - pressure[lev_idx+1]             ; pascal
       tot_layer_air_dens = - delta_pres / (mass_p_molec_air * mean_gravity) * 1e-4 ; molecule / cm^3
       
       tot_col_air = tot_col_air + tot_layer_air_dens

       ;; Calculate species mean layer density      
       species_mr_m = profile[lev_idx] + (profile[lev_idx + 1] - profile[lev_idx]) * interp_f
       species_dens_m = species_mr_m * tot_layer_air_dens
       
       tot_col_species = tot_col_species + species_dens_m

   endfor

   Return, tot_col_species  ; molecule / cm^2

END
