PRO Convert_ODell_Radiance, save_file, output_dir, frame_id_list, average=average_spec, l1a_file=l1a_hdf_filename, snr_coefs_file=snr_coefs_filename

   ; Make sure these are Long or else
   ; will get errors when dealing w/ large sav files
   average_len = 237L ; = (79 secs for FTS integration) / (1/3 sec for OCO integration)
   overlap_len = 15L  ; about 15 seconds between FTS measurements

   override_vals = { BadPixelMapVersionNum:  '-32767 -32767 -32767', $
                     sounding_latitude:      '34.2023', $
                     sounding_longitude:     '-118.1740', $
                     sounding_altitude:      '379', $
                     sounding_zenith:        '0.0', $
                     sounding_azimuth:       '0.0', $
                     spacecraft_alt:         '0.0', $
                     spacecraft_lat:         '0.0', $
                     spacecraft_lon:         '0.0' $
                   }
   override_tags = tag_names(override_vals)

   ; Load save file into memory
   print, 'Loading: ', save_file
   restore, save_file

   ; Load L1A file
   if n_elements(sound) ne 0 then begin
       print, 'Using internal sounding info: sound'
       sounding_ids      = sound.sounding_id
       frame_ids         = sound.frame_id
       frame_time_stamps = sound.frame_time_stamp

       ;sounding_latitude      = sound.sounding_latitude
       ;sounding_longitude     = sound.sounding_longitude
       ;sounding_altitude      = sound.sounding_altitude
       ;sounding_solar_zenith  = sound.sounding_solar_zenith
       ;sounding_solar_azimuth = sound.sounding_solar_azimuth
   endif else if n_elements(oco) ne 0 then begin
       print, 'Using internal sounding info: oco.sound'
       sounding_ids      = oco.sound.sounding_id
       frame_ids         = oco.sound.frame_id
       frame_time_stamps = oco.sound.frame_time_stamp

       ;sounding_latitude      = oco.sound.sounding_latitude
       ;sounding_longitude     = oco.sound.sounding_longitude
       ;sounding_altitude      = oco.sound.sounding_altitude
       ;sounding_solar_zenith  = oco.sound.sounding_solar_zenith
       ;sounding_solar_azimuth = oco.sound.sounding_solar_azimuth
   endif else begin
       if n_elements(l1a_hdf_filename) eq 0 then begin
           Message, /CONTINUE, 'No internal sounding info and no l1a file specified to retrieve sounding info'
           stop
       endif

       print, 'No internal sounding information, loading l1a: ', l1a_hdf_filename
       l1a_obj = obj_new('oco_hdf_input_file', l1a_hdf_filename)
       
       sounding_ids      = l1a_obj->read('sounding_id')
       frame_ids         = l1a_obj->read('frame_id')
       frame_time_stamps = l1a_obj->read('frame_time_stamp')
   endelse

   ; Convert to double to avoid overflows
   print, 'Converting sav data to doubles'
   if n_elements(o2a) ne 0 then begin
       o2a  = double(o2a)
       wco2 = double(wco2)
       sco2 = double(sco2)
   endif else if n_elements(oco) ne 0 then begin
       o2a  = double(oco.o2a)
       wco2 = double(oco.wco2)
       sco2 = double(oco.sco2)
   endif else begin
       Message, 'Could not find radiance data'
   endelse

   help, o2a, wco2, sco2

   ;sounding_tags = tag_names(sounding[0])
   ;print, 'Sounding Tags: ', sounding_tags

   num_frames = n_elements(frame_ids)
   num_footprints = 8
   num_pixels = (size(o2a))[1]

   if n_elements(snr_coefs_filename) gt 0 then $
     restore, snr_coefs_filename
      
   if n_elements(o2a_uncertainty) eq 0 then begin
       if n_elements(snr_coefs) eq 0 then begin
           Message, /CONTINUE, 'No snr coefs files supplied so can not calculate uncertainty'
           stop
       endif

       print, 'Calculating O2A uncertainty'
       o2a_uncertainty = calculate_oco_radiance_uncertainty(o2a, snr_coefs, 0)
   endif else begin
       o2a_uncertainty  = double(o2a_uncertainty)
   endelse

   if n_elements(wco2_uncertainty) eq 0 then begin
       print, 'Calculating WCO2 uncertainty'
       wco2_uncertainty = calculate_oco_radiance_uncertainty(wco2, snr_coefs, 1)
   endif else begin
       wco2_uncertainty = double(wco2_uncertainty)
   endelse

   if n_elements(sco2_uncertainty) eq 0 then begin
       print, 'Calculating SCO2 uncertainty'
       sco2_uncertainty = calculate_oco_radiance_uncertainty(sco2, snr_coefs, 2)
   endif else begin
       sco2_uncertainty = double(sco2_uncertainty)
   endelse

   help, o2a_uncertainty, wco2_uncertainty, sco2_uncertainty

   last_end_idx = -999
   for frame_idx = 0, num_frames-1 do begin
       frame_id = strtrim(string(frame_ids[frame_idx]),2)

       where_id = where(frame_ids eq frame_id)
       if(where_id[0] ge 0) then begin

           beg_frame_idx = frame_idx
           if keyword_set(average_spec) then begin
               end_frame_idx = beg_frame_idx + average_len

               if last_end_idx ge (beg_frame_idx - overlap_len) then $
                 Message, 'Frame ID: ' + frame_id + ' overlaps with frame ids used in last averaging'

               if end_frame_idx ge num_frames then $
                 Message, 'Frame ID: ' + frame_id + ' can not be averaged since it exceeds the bounds of available data'

               last_end_idx = end_frame_idx
           endif else begin
               end_frame_idx = frame_idx
           endelse
          
           ;; Use frame_id for time_stamp until it is fixed
           frame_time_stamp = strtrim(string(frame_time_stamps[frame_idx]),2)

          
           
           ;; Calculate azimuth/zenith for date/time
           if n_elements(sounding_solar_zenith) eq 0 or n_elements(sounding_solar_azimuth) eq 0 then begin
               zenith=0.0D0
               azimuth=0.0D0              
               solfac=0.0D0

               year   = fix(strmid(frame_time_stamp,  0, 4))
               month  = fix(strmid(frame_time_stamp,  5, 2))
               day    = fix(strmid(frame_time_stamp,  8, 2))
               hour   = fix(strmid(frame_time_stamp, 11, 2))
               minute = fix(strmid(frame_time_stamp, 14, 2))
               sec    = double(strmid(frame_time_stamp, 17, 6))

               ;; Calculate julian day of year
               doy=julday(month, day, year)-julday(1, 1, year)+1
           
               zensun, doy, hour+minute/60.0+sec/3600.0, $
                       override_vals.sounding_latitude, $
                       override_vals.sounding_longitude, $
                       zenith, azimuth, solfac
           endif

           for foot_idx = 0, num_footprints-1 do begin
               out_matrix = obj_new('matrix_file')

               sounding_id = sounding_ids[foot_idx, frame_idx]

               ;; Override with static values
               for key_idx = 0, n_elements(override_tags)-1 do begin
                   out_matrix->Set_Header_Keyword, strlowcase(override_tags[key_idx]), strtrim(string(override_vals.(key_idx)),2)
               endfor

               ;; Set these if they exist, replacing any overrides
               if n_elements(sounding_latitude) ne 0 then $
                 out_matrix->Set_Header_Keyword, 'sounding_latitude',  sounding_latitude[foot_idx, frame_idx]

               if n_elements(sounding_longitude) ne 0 then $
                 out_matrix->Set_Header_Keyword, 'sounding_longitude', sounding_longitude[foot_idx, frame_idx]

               if n_elements(sounding_altitude) ne 0 then $
                 out_matrix->Set_Header_Keyword, 'sounding_altitude',  sounding_altitude[foot_idx, frame_idx]

               if n_elements(sounding_solar_zenith) ne 0 then $
                 out_matrix->Set_Header_Keyword, 'sounding_solar_zenith', sounding_solar_zenith[foot_idx, frame_idx] $
               else $
                 out_matrix->Set_Header_Keyword, 'sounding_solar_zenith', strtrim(string(zenith),2)

               if n_elements(sounding_solar_azimuth) ne 0 then $
                 out_matrix->Set_Header_Keyword, 'sounding_solar_azimuth', sounding_solar_azimuth[foot_idx, frame_idx] $
               else $
                 out_matrix->Set_Header_Keyword, 'sounding_solar_azimuth', strtrim(string(azimuth),2)
               
               out_matrix->Set_Header_Keyword, 'sounding_id', strtrim(string(sounding_id),2)
               out_matrix->Set_Header_Keyword, 'frame_id', strtrim(string(frame_id),2)
               out_matrix->Set_Header_Keyword, 'frame_time_stamp', strtrim(string(frame_time_stamp),2)
               

               ;; Append data into a single matrix
               new_data = dblarr(2, num_pixels*3)

               
               ;; Sum pixels for all frames being averaged over
               for avg_frame_idx = beg_frame_idx, end_frame_idx do begin
                   new_data[0, 0:num_pixels-1] = new_data[0, 0:num_pixels-1] + $
                     o2a[*, foot_idx, avg_frame_idx]
                   new_data[0, num_pixels:(2*num_pixels)-1] = new_data[0, num_pixels:(2*num_pixels)-1] + $
                     wco2[*, foot_idx, avg_frame_idx]
                   new_data[0, (2*num_pixels):(3*num_pixels)-1] = new_data[0, (2*num_pixels):(3*num_pixels)-1] + $
                     sco2[*, foot_idx, avg_frame_idx]
                   
                   new_data[1, 0:num_pixels-1] = new_data[1, 0:num_pixels-1] + $
                     o2a_uncertainty[*, foot_idx, avg_frame_idx]^2
                   new_data[1, num_pixels:(2*num_pixels)-1] = new_data[1, num_pixels:(2*num_pixels)-1] + $
                     wco2_uncertainty[*, foot_idx, avg_frame_idx]^2
                   new_data[1, (2*num_pixels):(3*num_pixels)-1] = new_data[1, (2*num_pixels):(3*num_pixels)-1] + $
                     sco2_uncertainty[*, foot_idx, avg_frame_idx]^2
               endfor

               ; Compute average or error noise value for each pixel
               for pix_idx = 0, (num_pixels*3)-1 do begin
                   new_data[0, pix_idx] = new_data[0, pix_idx] / (end_frame_idx - beg_frame_idx + 1)
                   new_data[1, pix_idx] = sqrt(new_data[1, pix_idx])
               endfor

               out_matrix->set_data, new_data

               ;; Create filename based on sounding id               
               out_filename = string(output_dir, strtrim(string(sounding_id),2), FORMAT='(A0, "/", A0, ".dat")')

               print, 'Writing: ', out_filename
               out_matrix->write, out_filename, COLUMN_FORMAT=Replicate('E15.7', 2)
               ;; Release allocated memory
               new_data = 0
               out_matrix = 0
           endfor ;; foot_idx

       endif ;; if frame exists
   endfor ;; frame id loop

   ;; Just in case memory cleanup
   o2a = 0
   wco2 = 0
   sco2 = 0
   sounding = 0
   o2a_uncertainty = 0
   wco2_uncertainty = 0
   sco2_uncertainty = 0
END
