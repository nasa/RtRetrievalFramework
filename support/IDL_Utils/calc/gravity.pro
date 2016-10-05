function gravity, gdlat, altit

; Computes the effective Earth gravity at a given latitude and altitude.
; This is the sum of the gravitational and centripital accelerations.
; These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
; The Earth is assumed to be an oblate ellipsoid, with a ratio of the
; major to minor axes = sqrt(1+con) where con=.006738
; This eccentricity makes the Earth's gravititational field smaller at the
; poles and larger at the equator than if the Earth were a sphere of the
; same mass. It also makes the local mid-latitude gravity field not point
; toward the center of mass.
;
; Input Parameters:
;     gdlat       r*4  GeoDetric Latitude (degrees)
;     altit       r*4  Geometric Altitude (km)
;
; Output Parameter:
;     gravity     r*4  Effective Gravitational Acceleration (m/s2)
;
; Interestingly, since the centripital effect of the Earth's rotation
; (-ve at equator, 0 at poles) has almost the opposite shape to the
; second order gravitational field (+ve at equator, -ve at poles), their
; sum is almost constant so that the surface gravity can be approximated
; (to .07%) by the simple expression g = 0.99746*GM/radius**2, the latitude
; variation coming entirely from the variation of surface r with latitude.

   d2r   = 0.0D0
   gm    = 0.0D0
   omega = 0.0D0
   con   = 0.0D0
   shc   = 0.0D0
   eqrad = 0.0D0
   
   gclat  = 0.0D0 ; geocentric latitude (radians)
   radius = 0.0D0 ; radial distance (metres)
   
   ; scratch variables
   ff  = 0.0D0
   hh  = 0.0D0
   gee = 0.0D0

   d2r=3.14159265/180 ; Conversion from degrees to radians
   gm=3.9862216e+14   ; Gravitational constant times Earth's Mass (m3/s2)
   omega=7.292116E-05 ; Earth's angular rotational velocity (radians/s)
   con=.006738        ; (a/b)^2-1 where a & b are equatorial & polar radii
   shc=1.6235e-03     ; 2nd harmonic coefficient of Earth's gravity field 
   eqrad=6378178.     ; Equatorial Radius (m)

   ;; Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
   gclat=atan(tan(d2r*gdlat)/(1+con))  ; radians

   ;; On computers which crash at the poles try the following expression
   ;; gclat=d2r*gdlat-con*sin(d2r*gdlat)*cos(d2r*gdlat)/(1+con*cos(d2r*gdlat)^2)
   radius=1000*altit+eqrad/sqrt(1.+con*sin(gclat)^2)
   ff=(radius/eqrad)^2
   hh=radius*omega^2
   gee=gm/eqrad^2                       ; = gravity at Re
   gravity=(gee*(1-shc*(3*sin(gclat)^2-1)/ff)/ff-hh*cos(gclat)^2) $
           *(1+0.5*(sin(gclat)*cos(gclat)*(hh/gee+2*shc/ff^2))^2)

   return, gravity
end
