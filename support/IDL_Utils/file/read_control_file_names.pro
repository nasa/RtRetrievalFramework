FUNCTION Read_Control_File_Names, control_file, SV_SIZES=sv_sizes, GROUP_AEROSOL=group_aerosol

   ;; Rename certain items to a more appropriate name
   sv_name_replace   = [ [ 'PRESS',   'PSURF' ] ]
   use_element_name  = [ 'CO2', 'H2O', 'O2', 'ALBEDO', 'WINDSPEED' ]
   
   aerosol_element_names = [ 'CONT', 'OCEAN', 'IC_', 'WC_', 'LOW.?_', 'MID.?_', 'HIGH.?_' ]

   if not keyword_set(group_aerosol) then $
     use_element_name = [use_element_name, aerosol_element_names]

   ;; Read control_file
   OpenR, lun1, control_file, /GET_LUN

   parse_flag_dims   = 0
   end_sec_flag_dims = 0

   parse_flag_elem   = 0
   end_sec_flag_elem = 0
   need_elems        = 1

   line_str = ''

   ;; Grow as needed
   num_sv_arr   = 0
   sv_sizes     = intarr(1)
   sv_names_all = strarr(1, 500)
   sv_names_tmp = strarr(500)

   elem_idx = 0

   sv_name_idx = 0
   while not eof(lun1) do begin
       Readf, lun1, line_str

       if (strpos(line_str, 'SV Dimensions') ge 0) then begin
           parse_flag_dims = 1
       endif else if (strpos(line_str, 'SV_length') ge 0) then begin
           end_sec_flag_dims = 1
           parse_flag_dims = 0
       endif else if (parse_flag_dims) then begin
           line_parts = strsplit(line_str, ':', /EXTRACT)

           curr_name = strupcase( strcompress(strtrim(line_parts[0], 2)) )
           num_elem  = fix( strcompress(line_parts[1]) )

           space_pos = strpos(curr_name, ' ')
           while space_pos ne -1 do begin
               strput, curr_name, '_', space_pos
               space_pos = strpos(curr_name, ' ')
           endwhile
           
           ; Remove space
           replace_idx = where(strupcase(sv_name_replace[0, *]) eq curr_name )
           if replace_idx[0] ne -1 then $
             curr_name = sv_name_replace[1, replace_idx]

           ; Grow temp names array if needed
           if n_elements(sv_names_tmp) lt sv_name_idx + num_elem then begin
               sv_name_new = strarr( n_elements(sv_names_tmp) + num_elem )
               sv_name_new[0:n_elements(sv_names_tmp)-1] = sv_names_tmp
               sv_names_tmp = sv_name_new
           endif

           for sv_set_idx = sv_name_idx, sv_name_idx + num_elem - 1 do begin
               sv_names_tmp[sv_set_idx] = curr_name
           endfor

           sv_name_idx = sv_name_idx + num_elem
       endif else if (end_sec_flag_dims) then begin
           end_sec_flag_dims = 0
           
           where_sv_size = where(sv_sizes eq sv_name_idx)
           if where_sv_size[0] eq -1 then begin

               ; Increase number of sv arrangements
               num_sv_arr = num_sv_arr + 1

               ; grow sv_all, sv_names if needed
               num_sv_rows = n_elements( sv_names_all(0, *) )

               if n_elements(sv_sizes) lt num_sv_arr then begin
                   ;; should only happen num_sv_arr >= 2
                   sv_sizes_new = intarr(num_sv_arr)

                   sv_sizes_new(0:num_sv_arr-2) = sv_sizes(*)
                   sv_sizes = sv_sizes_new

                   num_sv_rows_new = num_sv_rows
                   if num_sv_rows lt sv_name_idx then $
                     num_sv_rows_new = sv_name_idx
                     
                   sv_names_new = strarr( num_sv_arr, num_sv_rows_new )
                   for col_idx = 0, num_sv_arr-2 do begin
                       sv_names_new(col_idx, 0:num_sv_rows-1) = sv_names_all(col_idx, 0:num_sv_rows-1)
                   endfor
                   sv_names_all = sv_names_new
               endif

               ; store data 
               sv_sizes[ num_sv_arr - 1 ] = sv_name_idx
               sv_names_all[ num_sv_arr - 1, 0:sv_name_idx-1] = sv_names_tmp[0:sv_name_idx-1]

               ; element names instead of dim specifiers
               sv_elems = strarr(sv_name_idx)
           endif 

           sv_name_idx = 0
       endif
       
       if (need_elems and strpos(line_str, 'SV Elements') ge 0) then begin
           parse_flag_elem = 1
       endif else if (parse_flag_elem and strpos(line_str, '#') ge 0) then begin
           end_sec_flag_elem = 1
           parse_flag_elem = 0
       endif else  if (parse_flag_elem) then begin
           line_parts = strsplit(line_str, /EXTRACT)
           
           curr_name = strupcase( strcompress(strtrim(line_parts[0], 2)) )
           sv_elems[elem_idx] = curr_name

           elem_idx = elem_idx + 1

       endif else if (end_sec_flag_elem) then begin
           end_sec_flag_elem = 0
           
           sv_names_uniq = sv_names_all[num_sv_arr-1, uniq(sv_names_all[num_sv_arr-1, 0:sv_sizes[num_sv_arr-1]-1])]
           sv_elems_uniq = sv_elems[uniq(sv_elems)]

           need_elems = 0
       endif
              
   endwhile
   Free_Lun, lun1

   max_sv_size = max(sv_sizes)
   sv_names_new = strarr(num_sv_arr, max_sv_size)
   for col_idx = 0, num_sv_arr-1 do begin
       sv_names_new[col_idx, 0:max_sv_size-1] = sv_names_all[col_idx, 0:max_sv_size-1]

       ;; Correct elements in sv_names_all using sv_elems
       for sv_idx = 0, sv_sizes[col_idx]-1 do begin
           where_elem_idx = where(sv_names_all[col_idx, sv_idx] eq sv_names_uniq)

           for chk_idx = 0, n_elements(use_element_name)-1 do begin
               if stregex(sv_elems_uniq[where_elem_idx[0]], '^'+use_element_name[chk_idx], /BOOLEAN, /FOLD_CASE) then begin
                   sv_names_new[col_idx, sv_idx] = sv_elems_uniq[where_elem_idx]
                   break
               endif
           endfor
       endfor

   endfor
   sv_names_all = sv_names_new
               
   Return, sv_names_all

END
