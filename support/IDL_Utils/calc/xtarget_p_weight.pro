FUNCTION XTarget_P_Weight, pressure, profile, PSURF=psurf

   x_target = 0.0E0

   n_levels = n_elements(profile)

   ;; Interpolate mixing ratio to surface pressure
   if n_elements(psurf) ne 0 then begin
       profile(n_levels-1) =                                           $
         (profile(n_levels-1) - profile(n_levels-2))                   $
         / (alog(pressure(n_levels-1))                              $
            - alog(pressure(n_levels-2)))                            $
         * (alog(psurf) - alog(pressure(n_levels-2)))  $
         + profile(n_levels-2)
       pressure(n_levels-1) = psurf
   endif

   ;; Pressure weighting function
   totp = pressure(n_levels-1) - pressure(0)
   for lev_idx = 0, n_levels-2 do begin
       dpres = pressure(lev_idx + 1) - pressure(lev_idx)
       press_wf_lev = dpres / totp
       
       rp = pressure(lev_idx + 1) / pressure(lev_idx)
       interp_f = pressure(lev_idx + 1) / dpres - 1.0D / alog(rp)
       
       target_mr_m = profile(lev_idx) + (profile(lev_idx + 1) - profile(lev_idx)) * interp_f
       
       x_target = x_target + press_wf_lev * target_mr_m
   endfor

   Return, x_target

END
