function compute_xco2, co2, p, t, h2o, ak=ak

; COMPUTES XCO2 from Profiles of CO2 and P on levels [TOA->surface ordering]
; If T, H2O are present, the water vapor correction will be made.
; (gravity must be computed in this case)
;
; Can take into account optional averaging kernel (the result will be unnormalized!)
; AK is the averaging kernel on levels

; co2 : profile of co2 on levels
; p   : pressure profile [Pa] on levels
; t   : temp profile [K] on levels
; h2o : h2o volume concentration profile on levels

nlev = n_elements(co2)
; Compute h2o correction if desired
if n_elements(t) eq nlev AND n_elements(h2o) eq nlev then begin
    epsilon = 0.622
    q = epsilon * h2o/ (1. + epsilon * h2o)
    qbar = (q[1:*] + q)*0.5
    gbar = gravity_profile(p,T,q,0., 45.)
    corr = (1-qbar)/gbar
endif else begin
    corr = fltarr(nlev-1) + 1.
endelse

; Get pressure weighting function on layers
dp = P[1:*]-P
wgt = dp * corr
wgt = wgt / total(wgt)

; Get the pressure weighting function on levels
wgtlev = fltarr(nlev)
wgtlev[0] = 0.5 * wgt[0]
wgtlev[1:nlev-2] = 0.5 * (wgt[1:nlev-2] + wgt[0:nlev-3] )
wgtlev[nlev-1] = 0.5 * wgt[nlev-2]

; compute XCO2
if n_elements(ak) eq 0 then ak = fltarr(nlev) + 1.
xco2 = total(wgtlev * co2 * ak)
return, xco2

END



