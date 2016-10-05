pro zensun,day,time,lat,lon,zenith,azimuth,solfac,sunrise=sunrise, $
           sunset=sunset,local=local,latsun=latsun,lonsun=lonsun
;+
; ROUTINE:      zensun
;
; PURPOSE:      Compute solar position information as a function of
;               geographic coordinates, date and time.
;               
; USEAGE:       zensun,day,time,lat,lon,zenith,azimuth,solfac,sunrise,sunset,
;                  local=local
;
; INPUT:
;   day         Julian day (positive scalar or vector)
;               (spring equinox =  80)
;               (summer solstice= 171)
;               (fall equinox   = 266)
;               (winter solstice= 356)
;
;   time        Universal Time in hours (scalar or vector) 
;
;   lat         geographic latitude of point on earth's surface (degrees)
;
;   lon         geographic longitude of point on earth's surface (degrees)
;
; OUTPUT:
;
;   zenith      solar zenith angle (degrees)
;
;   azimuth     solar azimuth  (degrees) 
;               Azimuth is measured clockwise from due north 
;
;   solfac      Solar flux multiplier.  SOLFAC=cosine(ZENITH)/RSUN^2
;               where rsun is the current earth-sun distance in
;               astronomical units.
;
;               NOTE: SOLFAC is negative when the sun is below the horizon 
;
;
; KEYWORD INPUT:
;
;   local       if set, TIME is specified as a local time and SUNRISE
;               and SUNSET are output in terms of local time
; 
;               NOTE: "local time" is defined as UT + local_offset
;                     where local_offset is fix((lon+sign(lon)*7.5)/15)
;                     with -180 < lon < 180
;
;                     Be aware, there are no fancy provisions for
;                     day-light savings time and other such nonsense.
;
; KEYWORD OUTPUT:
;
;   sunrise     Time of sunrise (hours)
;   sunset      Time of sunset  (hours)
;
;   latsun      the latitude of the sub-solar point (fairly constant over day)
;               Note that daily_minimum_zenith=abs(latsun-lat)
;
;   lonsun      the longitude of the sub-solar point
;               Note that at 12 noon local time (lon-lonsun)/15. is the
;               number of minutes by which the sun leads the clock.
;
;;              Often used           lat,   lon 
;
;               Santa Barbara:        34.410,-119.727
;               SGP Cart Site:        36.605,-97.485
;               North Slope:          69.318,-151.015
;               Palmer Station:       -64.767,-64.067
; 
;; EXAMPLE 1:   Compute the solar flux at Palmer Station for day 283
;
;               time=findgen(1+24*60)/60
;               zensun,283,time,-64.767,-64.067,z,a,sf
;               solflx=sf*s
;               plot,time,solflx
;
;               where s is the TOA solar irradiance at normal incidence:
;
;               s=1618.8   ; W/m2/micron for AVHRR1 GTR100 
;               s= 976.9   ; W/m2/micron for AVHRR2 GTR100
;               s=1685.1   ; W/m2/micron for 410nm GTR100
;               s= 826.0   ; W/m2/micron for 936nm GTR100
;               s=1.257e17 ; photons/cm2/s PAR GTR100
;               s=1372.9   ; w/m2 total broadband
;
;
;; EXAMPLE 2:   Find time of sunrise and sunset for current day
; 
;     doy=julday()-julday(1,0,1994)                          
;     zensun,doy,12,34.456,-119.813,z,a,s,sunrise=sr,sunset=ss,/local &$
;     zensun,doy,[sr,.5*(sr+ss),ss],34.456,-119.813,z,az,/local &$
;     print,'sunrise: '+hms(3600*sr)+      ' PST   azimuth angle: ',az(0)  &$
;     print,'sunset:  '+hms(3600*ss)+      ' PST   azimuth angle: ',az(2)  &$
;     print,'zenith:  '+hms(1800*(ss+sr))+ ' PST   zenith angle:  ',z(1) 
;   
;  AUTHOR:      Paul Ricchiazzi        23oct92
;               Earth Space Research Group,  UCSB
; 
;  REVISIONS:
;  
; jan94: use spline fit to interpolate on EQT and DEC tables
; jan94: output SUNRISE and SUNSET, allow input/output in terms of local time 
; jan97: declare eqtime and decang as floats.  previous version 
;         this caused small offsets in the time of minimum solar zenith
;-
; PROCEDURE: 
;
; 1.  Calculate the subsolar point latitude and longitude, based on
;     DAY and TIME. Since each year is 365.25 days long the exact
;     value of the declination angle changes from year to year.  For
;     precise values consult THE AMERICAN EPHEMERIS AND NAUTICAL
;     ALMANAC published yearly by the U.S. govt. printing office.  The
;     subsolar coordinates used in this code were provided by a
;     program written by Jeff Dozier.
;
;  2. Given the subsolar latitude and longitude, spherical geometry is
;     used to find the solar zenith, azimuth and flux multiplier.
;
;  eqt = equation of time (minutes)  ; solar longitude correction = -15*eqt
;  dec = declination angle (degrees) = solar latitude 
;
; LOWTRAN v7 data (25 points)
;     The LOWTRAN solar position data is characterized by only 25 points.
;     This should predict the subsolar angles within one degree.  For
;     increased accuracy add more data points.
;
;nday=[   1.,    9.,   21.,   32.,   44.,   60.,  91.,  121.,  141.,  152.,$
;       160.,  172.,  182.,  190.,  202.,  213., 244.,  274.,  305.,  309.,$
;       325.,  335.,  343.,  355.,  366.]
;
;eqt=[ -3.23, -6.83,-11.17,-13.57,-14.33,-12.63, -4.2,  2.83,  3.57,  2.45,$
;       1.10, -1.42, -3.52, -4.93, -6.25, -6.28,-0.25, 10.02, 16.35, 16.38,$
;       14.3, 11.27,  8.02,  2.32, -3.23]
;
;dec=[-23.07,-22.22,-20.08,-17.32,-13.62, -7.88, 4.23, 14.83, 20.03, 21.95,$
;      22.87, 23.45, 23.17, 22.47, 20.63, 18.23, 8.58, -2.88,-14.18,-15.45,$
;     -19.75,-21.68,-22.75,-23.43,-23.07]
;
; Analemma information from Jeff Dozier
;     This data is characterized by 74 points
;

if n_params() eq 0 then begin
  xhelp,'zensun'
  return
endif  

nday=[  1.0,   6.0,  11.0,  16.0,  21.0,  26.0,  31.0,  36.0,  41.0,  46.0,$
       51.0,  56.0,  61.0,  66.0,  71.0,  76.0,  81.0,  86.0,  91.0,  96.0,$
      101.0, 106.0, 111.0, 116.0, 121.0, 126.0, 131.0, 136.0, 141.0, 146.0,$
      151.0, 156.0, 161.0, 166.0, 171.0, 176.0, 181.0, 186.0, 191.0, 196.0,$
      201.0, 206.0, 211.0, 216.0, 221.0, 226.0, 231.0, 236.0, 241.0, 246.0,$
      251.0, 256.0, 261.0, 266.0, 271.0, 276.0, 281.0, 286.0, 291.0, 296.0,$
      301.0, 306.0, 311.0, 316.0, 321.0, 326.0, 331.0, 336.0, 341.0, 346.0,$
      351.0, 356.0, 361.0, 366.0]

eqt=[ -3.23, -5.49, -7.60, -9.48,-11.09,-12.39,-13.34,-13.95,-14.23,-14.19,$
     -13.85,-13.22,-12.35,-11.26,-10.01, -8.64, -7.18, -5.67, -4.16, -2.69,$
      -1.29, -0.02,  1.10,  2.05,  2.80,  3.33,  3.63,  3.68,  3.49,  3.09,$
       2.48,  1.71,  0.79, -0.24, -1.33, -2.41, -3.45, -4.39, -5.20, -5.84,$
      -6.28, -6.49, -6.44, -6.15, -5.60, -4.82, -3.81, -2.60, -1.19,  0.36,$
       2.03,  3.76,  5.54,  7.31,  9.04, 10.69, 12.20, 13.53, 14.65, 15.52,$
      16.12, 16.41, 16.36, 15.95, 15.19, 14.09, 12.67, 10.93,  8.93,  6.70,$
       4.32,  1.86, -0.62, -3.23]

dec=[-23.06,-22.57,-21.91,-21.06,-20.05,-18.88,-17.57,-16.13,-14.57,-12.91,$
     -11.16, -9.34, -7.46, -5.54, -3.59, -1.62,  0.36,  2.33,  4.28,  6.19,$
       8.06,  9.88, 11.62, 13.29, 14.87, 16.34, 17.70, 18.94, 20.04, 21.00,$
      21.81, 22.47, 22.95, 23.28, 23.43, 23.40, 23.21, 22.85, 22.32, 21.63,$
      20.79, 19.80, 18.67, 17.42, 16.05, 14.57, 13.00, 11.33,  9.60,  7.80,$
       5.95,  4.06,  2.13,  0.19, -1.75, -3.69, -5.62, -7.51, -9.36,-11.16,$
     -12.88,-14.53,-16.07,-17.50,-18.81,-19.98,-20.99,-21.85,-22.52,-23.02,$
     -23.33,-23.44,-23.35,-23.06]

;
; compute the subsolar coordinates
;

tt=((fix(day)+time/24.-1.) mod 365.25) +1.  ;; fractional day number
                                            ;; with 12am 1jan = 1.

if n_elements(tt) gt 1 then begin
  eqtime=tt-tt                              ;; this used to be day-day, caused 
  decang=eqtime                             ;; error in eqtime & decang when a
  ii=sort(tt)                               ;; single integer day was input
  eqtime(ii)=spline(nday,eqt,tt(ii))/60.    
  decang(ii)=spline(nday,dec,tt(ii))
endif else begin
  eqtime=spline(nday,eqt,tt)/60.
  decang=spline(nday,dec,tt)
endelse  
latsun=decang

if keyword_set(local) then begin
  lonorm=((lon + 360 + 180 ) mod 360 ) - 180.
  tzone=fix((lonorm+7.5)/15)
  index = where(lonorm lt 0, cnt)
  if (cnt gt 0) then tzone(index) = fix((lonorm(index)-7.5)/15)
  ut=(time-tzone+24.) mod 24.                  ; universal time
  noon=tzone+12.-lonorm/15.                    ; local time of noon
endif else begin
  ut=time
  noon=12.-lon/15.                             ; universal time of noon
endelse

lonsun=-15.*(ut-12.+eqtime)

; compute the solar zenith, azimuth and flux multiplier

t0=(90.-lat)*!dtor                            ; colatitude of point
t1=(90.-latsun)*!dtor                         ; colatitude of sun

p0=lon*!dtor                                  ; longitude of point
p1=lonsun*!dtor                               ; longitude of sun

zz=cos(t0)*cos(t1)+sin(t0)*sin(t1)*cos(p1-p0) ; up          \
xx=sin(t1)*sin(p1-p0)                         ; east-west    > rotated coor
yy=sin(t0)*cos(t1)-cos(t0)*sin(t1)*cos(p1-p0) ; north-south /

azimuth=atan(xx,yy)/!dtor                     ; solar azimuth 
zenith=acos(zz)/!dtor                         ; solar zenith

rsun=1.-0.01673*cos(.9856*(tt-2.)*!dtor)      ; earth-sun distance in AU
solfac=zz/rsun^2                              ; flux multiplier

if n_elements(time) eq 1 then begin
  angsun=6.96e10/(1.5e13*rsun)               ; solar disk half-angle
  ;angsun=0.  
  arg=-(sin(angsun)+cos(t0)*cos(t1))/(sin(t0)*sin(t1))
  sunrise = arg - arg 
  sunset  = arg - arg + 24.
  index = where(abs(arg) le 1, cnt)
  if (cnt gt 0) then begin
    dtime=acos(arg(index))/(!dtor*15)
    sunrise(index)=noon-dtime-eqtime(index)
    sunset(index)=noon+dtime-eqtime(index)
  endif
endif
return
end
