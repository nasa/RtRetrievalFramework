pro create_sd_info,runlog

; this routine takes a runlog file and create a soundinginfo file
; this soundinginfo file can then be used for the convolution of FTS spectra
; to obtain OCO files with the correct header information


lun1=1
lun2=2
str=''
openr,lun1,runlog
WHILE ~ EOF(1) DO BEGIN

  readf,lun1,str
  str1=strsplit(str,/extract)
  name=str1(0)

  name1=strsplit(name,/extract)

  month=strcompress(strmid(name1,6,2),/remove_all)
  day=strcompress(strmid(name1,8,2),/remove_all)

  year=strcompress(str1(1),/remove_all)

  dayofyear=long(str1(2))

  frac=float(str1(3))
  hours=long(frac)
  min=long((frac-float(hours))*60.)

  sec=long((frac-hours-float(min)/60.)*3600.)

  
  sec_frac=(frac-hours-float(min)/60.)*3600.-sec
  
  lat=str1(4)
  lon=str1(5)
  alt=str1(6)
  sza=str1(7)
  
 openw,lun2,'soundinginfo.dat_'+name
  
  timestamp=' frame_time_stamp = '+year+'-'+month+'-'+day+'T'
  h=string(hours,format='(I2.2)')
  m=string(min,format='(I2.2)')
  s=string(sec,format='(I2.2)')  

  s_frac=strcompress(string(sec_frac,format='(F5.3)'),/remove_all)
  s_fr=strmid(s_frac,1,4)
  time=h+':'+m+':'+s+s_fr+'Z'

  printf,lun2,timestamp+time
  printf,lun2,' sounding_altitude = ',alt 
  printf,lun2,' sounding_latitude = ',lat
  printf,lun2,' sounding_longitude = ',lon
  printf,lun2,' sounding_azimuth =  180.000000000000'
  printf,lun2,' sounding_zenith = ',sza
  close,lun2
endwhile

close,lun1




end
