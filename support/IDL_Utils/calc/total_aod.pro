function total_aod, pressure, aerosol, LAYER_AOD=layer_aod

  pre_dim  = size(pressure,  /DIMENSION)
  aer_dim  = size(aerosol,  /DIMENSION)

  ;; convert from log if aerosols are in log. 
  ;; If there are any negatives (should all since no aerosol > 1
  ;; then we have log retrieval of aerosol
  where_neg = where(aerosol lt 0.0, neg_count)
  if neg_count gt 0 then $
    aerosol = exp(aerosol)

  if n_elements(pre_dim) eq 1 then begin
      num_levels = min([aer_dim[1], pre_dim[0]])
  endif else begin
      num_levels = min([aer_dim[1], pre_dim[1]])
  endelse

  ;; compute total aod
  layer_aod = dblarr(aer_dim[1]-1)
  for aer_idx = 0, aer_dim[0] - 1 do begin
      for pres_idx = 0, num_levels - 2 do begin
          layer_aod[pres_idx] = layer_aod[pres_idx] + ((pressure[pres_idx+1] - pressure[pres_idx]) * $
                                                       ((aerosol[aer_idx, pres_idx] + aerosol[aer_idx, pres_idx + 1]) / 2.0D0))
      endfor
  endfor

  total_aod = total(layer_aod)

  return, total_aod

end
