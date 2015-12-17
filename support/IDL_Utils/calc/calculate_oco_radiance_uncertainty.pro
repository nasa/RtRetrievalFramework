function calculate_oco_radiance_uncertainty, radiance, snr_coefs, band_idx

   ;; Radiance : input radiance in photons/s/m^2/micron/sr
   ;;            dimension is [1016,8,nframes]
   ;; snr_coefs first two dimensions match dirst two dimensions of radiance

   max_ms = [1.4e21, 4.9e20, 1.7e20] ;; max meas signal by band

   NEN = radiance*0.0D0
   
   nframes = n_elements(radiance[0,0,*])
   
   for i = 0L, nframes-1 do begin
       NEN[*,*,i] = (max_ms[band_idx] / 100.0D0) * $
                    sqrt( abs(100.0D0 * radiance[*,*, i] / max_ms[band_idx]) * $
                          snr_coefs[*, *, band_idx, 0]^2 + snr_coefs[*, *, band_idx, 1]^2)
   endfor
   
   return, NEN

END
