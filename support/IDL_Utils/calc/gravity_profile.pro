function gravity_profile, P, T, q, elevation, lat

; p is ordered space to surface, in Pascals
; elevation is in km
; T is temperature [K] on levels
; q is specific humidity on levels
; lat is in degrees
;
; returns the gravity profile [m/s^2] on LAYERS

    efact = 0.60776532
    Rd = 287.04742
    nz = n_elements(P) - 1
    z = P*0.
    z[nz] = elevation
    grav = fltarr(nz)
    for i = nz-1, 0, -1 do begin
        qbar = (q[i] + q[i+1])*0.5
        Pbar = (p[i] + p[i+1])*0.5
        Tbar = (t[i] + t[i+1])*0.5
        zlower = z[i+1]
        dP = p[i+1] - p[i]
        Tv = Tbar * (1.0d0 + qbar * efact)
        logratio = alog(p[i+1] / p[i])
        grav[i] = jpl_gravity(lat, zlower)
        dz = logratio * Tv * Rd / grav[i] * 1d-3 ; km
        grav[i] = jpl_gravity(lat, zlower + 0.5*dz)
        dz = logratio * Tv * Rd / grav[i] * 1d-3 ; km
        z[i] = zlower + dz
    endfor

    return, grav

END
