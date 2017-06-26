;; given a key (from OCO simulated spectrum file), assign
;; corresponding L1B element
pro defineKeyIndex, table

  record = {translate, key:'', element:'', value:''}

  table = replicate({translate}, 61)

  table[ 1] = {key:'head_len_and_file_len', element:'UNKNOWN', value:''}
  table[ 2] = {key:'spectrum_name',         element:'UNKNOWN', value:''}
  table[ 3] = {key:'instrumentflag',        element:'UNKNOWN', value:''}
  table[ 4] = {key:'row_time',              element:'frame_time_stamp', value:''} ; 4.6.40
  table[ 5] = {key:'row_time_secs',         element:'frame_time', value:''} ; 4.6.39
  table[ 6] = {key:'row_inst_status',       element:'frame_inst_status', value:''} ; 4.6.37
  table[ 7] = {key:'row_err_status',        element:'frame_err_status', value:''} ; 4.6.36
  table[ 8] = {key:'row_qual_flag',         element:'frame_qual_flag', value:''} ; 4.6.38
  table[ 9] = {key:'x_pos',                 element:'x_pos', value:''} ; 4.6.130
  table[10] = {key:'y_pos',                 element:'y_pos', value:''} ; 4.6.133
  table[11] = {key:'z_pos',                 element:'z_pos', value:''} ; 4.6.135
  table[12] = {key:'x_vel',                 element:'x_vel', value:''} ; 4.6.131
  table[13] = {key:'y_vel',                 element:'y_vel', value:''} ; 4.6.134
  table[14] = {key:'z_vel',                 element:'z_vel', value:''} ; 4.6.136
  table[15] = {key:'roll',                  element:'roll', value:''} ; 4.6.94
  table[16] = {key:'pitch',                 element:'pitch', value:''} ; 4.6.66
  table[17] = {key:'yaw',                   element:'yaw', value:''} ; 4.6.132
  table[18] = {key:'spacecraft_lat',        element:'spacecraft_lat', value:''} ; 4.6.113
  table[19] = {key:'spacecraft_lon',        element:'spacecraft_lon', value:''} ; 4.6.114
  table[20] = {key:'spacecraft_alt',        element:'spacecraft_alt', value:''} ; 4.6.112
  table[21] = {key:'sounding_alt_geoid',    element:'footprint_alt_geoid', value:''} ; ???
  table[22] = {key:'sounding_lat_geoid',    element:'footprint_lat_geoid', value:''} ; 4.6.27
  table[23] = {key:'sounding_lon_geoid',    element:'footprint_lon_geoid', value:''} ; 4.6.29
  table[24] = {key:'sounding_alt_topo',     element:'footprint_alt_topo', value:''} ; ???
  table[25] = {key:'sounding_lat_topo',     element:'footprint_lat_topo', value:''} ; 4.6.26
  table[26] = {key:'sounding_lon_topo',     element:'footprint_lon_topo', value:''} ; 4.6.28
  table[27] = {key:'solar_azimuth_topo',    element:'solar_azimuth', value:''} ; 4.6.106
  table[28] = {key:'solar_zenith_topo',     element:'solar_zenith', value:''} ; 4.6.107
  table[29] = {key:'boresight_azimuth_topo',element:'satellite_azimuth', value:''} ; 4.6.96
  table[30] = {key:'boresight_zenith_topo', element:'satellite_zenith', value:''} ; 4.6.97
  table[31] = {key:'sounding_mode_flag',    element:'sounding_mode_flag', value:''}
  table[32] = {key:'sounding_qual_flag',    element:'sounding_qual_flag', value:''} ; 4.5.5
  table[33] = {key:'temp_instrument',       element:'temperature_internal', value:''} ; 4.5.5
  table[34] = {key:'press_instrument',      element:'pressure_internal', value:''} ; 4.5.5
  table[35] = {key:'hum_instrument',        element:'humidity_internal', value:''} ; 4.5.5
  table[36] = {key:'temp_outside',          element:'temperature_external', value:''} ; 4.5.5
  table[37] = {key:'press_outside',         element:'pressure_external', value:''} ; 4.5.5
  table[38] = {key:'hum_outside',           element:'humidity_external', value:''} ; 4.5.5
  table[39] = {key:'fov',                   element:'fov', value:''}
  table[40] = {key:'snr',                   element:'snr', value:''}
  table[41] = {key:'year',                  element:'year', value:''}
  table[42] = {key:'day',                   element:'day', value:''}
  table[43] = {key:'zpdtime',               element:'zpd_time_of_ut_day', value:''}
  table[44] = {key:'n_spec',                element:'n_spec', value:''}
  table[45] = {key:'n_zl_para',             element:'n_zl_para', value:''}
  table[46] = {key:'zero(1)',               element:'zero(1)', value:''}
  table[47] = {key:'zero(2)',               element:'zero(2)', value:''}
  table[48] = {key:'zero(3)',               element:'zero(3)', value:''}
  table[49] = {key:'n_ils_para',            element:'n_ils_para', value:''}
  table[50] = {key:'n_ils_wndepend',        element:'n_ils_wndepend', value:''}
  table[51] = {key:'ils(1)',                element:'ils(1)', value:''}
  table[52] = {key:'ils(2)',                element:'ils(2)', value:''}
  table[53] = {key:'ils(3)',                element:'ils(3)', value:''}
  table[54] = {key:'n_d_para',              element:'n_d_para', value:''}
  table[55] = {key:'n_pixel(1)',            element:'n_pixel(1)', value:''}
  table[56] = {key:'dispersion(1)',         element:'dispersion(1)', value:''}
  table[57] = {key:'n_pixel(2)',            element:'n_pixel(2)', value:''}
  table[58] = {key:'dispersion(2)',         element:'dispersion(2)', value:''}
  table[59] = {key:'n_pixel(3)',            element:'n_pixel(3)', value:''}
  table[60] = {key:'dispersion(3)',         element:'dispersion(3)', value:''}

end

;; Given an element in the table, return its index
function ElementIndex, table, element

  index = -1
  for i = 0, n_elements(table) - 1 do begin

      if (strcmp(strtrim(table[i].element), strtrim(element)) eq 1) then begin
          index = i
          break
      endif
      
  endfor

  if (index eq -1) then begin
      errMsg = "Can't find element " + element + " in table"
      message, errMsg
  endif

  return, index
end

;; Given a key in the table, return its index
function KeyIndex, table, key

  index = -1
  for i = 0, n_elements(table) - 1 do begin

      if (strcmp(strtrim(table[i].key), strtrim(key)) eq 1) then begin
          index = i
          break
      endif
      
  endfor

  if (index eq -1) then begin
      errMsg = "Can't find key " + key + " in table"
      message, errMsg
  endif

  return, index
end


;; Read the OCO simulated spectrum file (OCO_1_00001.0001).  Put the
;; header information in table[].value, and the radiances in arrays
pro readAsciiSpectra, filename, table, radiance_o2, radiance_o2_noise, $
                      radiance_weak_co2, radiance_weak_co2_noise,      $
                      radiance_strong_co2, radiance_strong_co2_noise,  $
                      npixels, status

  openr, lun, filename, error=status, /get_lun

  if (status ne 0) then begin
      print, "Can't open ", filename
      return
  endif

  key = ''
  value = ''
  for i = 2, n_elements(table) - 1 do begin
      readf, lun, key, value, format='(1x, a22, a39)'
      index = KeyIndex(table, key)
      table[index].value = value
;      print, i, ' ', table[index].element, table[index].value
  endfor

  npixels = intarr(3)

  index = KeyIndex(table, "n_pixel(1)")
  npixels[0] = table[index].value
  radiance_o2 = fltarr(npixels[0])
  radiance_o2_noise = fltarr(npixels[0])

  for i = 0, npixels[0] - 1 do begin
      readf, lun, rad, noise
      radiance_o2[i] = rad
      radiance_o2_noise[i] = noise
  endfor

  index = KeyIndex(table, "n_pixel(2)")
  npixels[1] = table[index].value
  radiance_weak_co2 = fltarr(npixels[1])
  radiance_weak_co2_noise = fltarr(npixels[1])

  for i = 0, npixels[1] - 1 do begin
      readf, lun, rad, noise
      radiance_weak_co2[i] = rad
      radiance_weak_co2_noise[i] = noise
  endfor

  index = KeyIndex(table, "n_pixel(3)")
  npixels[2] = table[index].value
  radiance_strong_co2 = fltarr(npixels[2])
  radiance_strong_co2_noise = fltarr(npixels[2])

  for i = 0, npixels[2] - 1 do begin
      readf, lun, rad, noise
      radiance_strong_co2[i] = rad
      radiance_strong_co2_noise[i] = noise
  endfor

  free_lun, lun

end


;; Write a string to the file, padded to 40 characters if needed
pro writeString, lun, string

  tail = '                                        ' ; 40 spaces

  stringToWrite = strmid(strtrim(string) + tail, 0, 40)
  printf, lun, stringToWrite

end

;; Add leading zeroes to the array to make it the specified size
pro padArray, array, new_size

  tmp = fltarr(new_size)
  bound2 = new_size - 1
  bound1 = bound2 - (n_elements(array) - 1)
  
  tmp[bound1:bound2] = array
  array = tmp

end

; writeBinarySpectra, description, ascii, binary
; for example,
;   writeBinarySpectra, 'al_Jan1_od01', 'spectra/ascii/al_Jan1_od01', 'spectra/binary/al_Jan1_od01.bin'
pro writeBinarySpectra, description, asciiFile, binaryFile

  LITTLE_ENDIAN = (byte (1, 0, 1))[0] ; true if this is a little endian machine

  ;; set up the key to element name mapping
  defineKeyIndex, table

  ;; read the simulated spectra file
  readAsciiSpectra, asciiFile, table, radiance_o2, radiance_o2_noise,  $
    radiance_weak_co2, radiance_weak_co2_noise, radiance_strong_co2,   $
    radiance_strong_co2_noise, npixels, status

  if (status ne 0) then return

  ;; add leading zeros to pad arrays to 1024 values
  padArray, radiance_o2, 1024
  padArray, radiance_o2_noise, 1024
  npixels[0] = 1024

  padArray, radiance_weak_co2, 1024
  padArray, radiance_weak_co2_noise, 1024
  npixels[1] = 1024

  padArray, radiance_strong_co2, 1024
  padArray, radiance_strong_co2_noise, 1024
  npixels[2] = 1024


  numColors = 1024L
  numDispersionCoeffs = 2L
  numFrames = 1L
  numILSCoeffs = 10L
  numSoundings = 1L
  numSpatialRows = 220L
  numSpectrometers = 3L

  openw, lun, binaryFile, /get_lun

;; File description
  writeString, lun, description

;; Name
  writeString, lun, 'Hari Nair'

;; File creation date
  writeString, lun, systime(/UTC) + ' UTC'

;; Machine description
  spawn, 'uname -niop', machine
  writeString, lun, machine

;; Begin Dimensions
  writeString, lun, 'Begin Dimensions'

  writeString, lun, 'Sounding'
  printf, lun, numSoundings

  writeString, lun, 'Frame'
  printf, lun, numFrames

  writeString, lun, 'FPAColor'
  printf, lun, numColors

;; End Dimensions
  writeString, lun, 'End Dimensions'

;; Begin Data
  writeString, lun, 'Begin Data'

  ignore = 0
  if (ignore ne 0) then begin
      element = 'footprint_alt_topo'
      writeString, lun, element
      writeString, lun, 'm'
      printf, lun, 4            ; bytes per element
      printf, lun, numSoundings * numFrames ; number of elements
      index = ElementIndex(table, element)
      writeu, lun, table[index].value
      printf, lun, ''           ; new line

      element = 'footprint_lat_topo'
      writeString, lun, element
      writeString, lun, 'deg'
      printf, lun, 4            ; bytes per element
      printf, lun, numSoundings * numFrames ; number of elements
      index = ElementIndex(table, element)
      writeu, lun, table[index].value
      printf, lun, ''           ; new line

      element = 'footprint_lon_topo'
      writeString, lun, element
      writeString, lun, 'deg'
      printf, lun, 4            ; bytes per element
      printf, lun, numSoundings * numFrames ; number of elements
      index = ElementIndex(table, element)
      writeu, lun, table[index].value
      printf, lun, ''           ; new line

      element = 'solar_azimuth'
      writeString, lun, element
      writeString, lun, 'deg'
      printf, lun, 4            ; bytes per element
      printf, lun, numSoundings * numFrames ; number of elements
      index = ElementIndex(table, element)
      writeu, lun, table[index].value
      printf, lun, ''           ; new line

      element = 'solar_zenith'
      writeString, lun, element
      writeString, lun, 'deg'
      printf, lun, 4            ; bytes per element
      printf, lun, numSoundings * numFrames ; number of elements
      index = ElementIndex(table, element)
      writeu, lun, table[index].value
      printf, lun, ''           ; new line

      ils = fltarr(numSpectrometers)

      key = 'ils(1)'
      index = KeyIndex(table, key)
      ils[0] = table[index].value

      key = 'ils(2)'
      index = KeyIndex(table, key)
      ils[1] = table[index].value

      key = 'ils(3)'
      index = KeyIndex(table, key)
      ils[2] = table[index].value

      element = 'instrument_line_shape_coeff'
      writeString, lun, element
      writeString, lun, ''      ; units?
      printf, lun, 4            ; bytes per element
      printf, lun, numSpectrometers * numSpatialRows * numColors * numILSCoeffs ; number of elements
      for j = 1L, numILSCoeffs do begin
          for i = 1L, numSpatialRows * numColors do begin
              writeu, lun, ils
          endfor
          ils = [0., 0., 0.] ; values from the ASCII file are coefficient 1, set the rest to 0
      endfor
      printf, lun, ''           ; new line

  endif

  dispersion = fltarr(numDispersionCoeffs * numSpectrometers) ; 2 * 3

  key = 'dispersion(1)'
  index = KeyIndex(table, key)
  dispString = table[index].value
  reads, dispString, c0, c1
  dispersion[0] = c0
  dispersion[1] = c1

  key = 'dispersion(2)'
  index = KeyIndex(table, key)
  dispString = table[index].value
  reads, dispString, c0, c1
  dispersion[2] = c0
  dispersion[3] = c1

  key = 'dispersion(3)'
  index = KeyIndex(table, key)
  dispString = table[index].value
  reads, dispString, c0, c1
  dispersion[4] = c0
  dispersion[5] = c1

  element = 'DispersionCoefSamp'
  writeString, lun, element
  writeString, lun, '' ; units?
  printf, lun, 4                ; bytes per element
  printf, lun, numSpectrometers * numDispersionCoeffs ; number of elements
  writeu, lun, dispersion
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_o2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[0] ; number of elements
  writeu, lun, radiance_o2
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_uncertainty_o2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[0] ; number of elements
  writeu, lun, radiance_o2_noise
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_weak_co2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[1] ; number of elements
  writeu, lun, radiance_weak_co2
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_uncertainty_weak_co2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[1] ; number of elements
  writeu, lun, radiance_weak_co2_noise
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_strong_co2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[2] ; number of elements
  writeu, lun, radiance_strong_co2
  printf, lun, ''               ; new line

  writeString, lun, 'radiance_uncertainty_strong_co2'
  writeString, lun, 'photons sec-1 m-2 sr-1 micrometer-1'
  printf, lun, 4                ; bytes per element
  printf, lun, numSoundings * numFrames * npixels[2] ; number of elements
  writeu, lun, radiance_strong_co2_noise
  printf, lun, ''               ; new line

;; End Data
  writeString, lun, 'End Data'

  free_lun, lun

end

pro writeAll

  for i = 1, 261 do begin
      dir = string(i, format='("sim_",i3.3)')
      ascii = 'spectra/one_orbit/ascii/' + dir + '/spec/OCO/OCO_1_00001.0001'
      caldat, systime(/UTC, /julian), m, d, y
      date = string(y-2000,m,d, format='("_", 3(I2.2))')
      binary = 'spectra/one_orbit/binary/' + dir + date + '.bin'

      writeBinarySpectra, dir, ascii, binary
  endfor

end
