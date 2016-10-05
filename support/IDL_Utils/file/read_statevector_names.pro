FUNCTION Read_Statevector_Names, sv_names_file, GROUP_AEROSOL=group_aerosol

   aerosol_element_names = [ 'CONT$', 'OCEAN$', 'CLAT_?', 'IC_?', 'WC_?', 'LOW.?$', 'MID.?$', 'HIGH.?$', 'AERO.?$' ]
   aerosol_group_format = '("AEROSOL_",I0)'

   name_obj = Obj_New('MATRIX_FILE', sv_names_file)

   names = (name_obj->get_data())[1,*]


   curr_aer_name = ''
   curr_aer_idx  = 0
   for name_idx = 0, n_elements(names)-1 do begin
       names[name_idx] = strupcase(names[name_idx])

       if stregex(names[name_idx], '_[0-9]{3}$', /BOOLEAN) then begin
           names[name_idx] = strmid(names[name_idx], 0, strlen(names[name_idx])-4)
       endif

       if Keyword_Set(group_aerosol) then begin
           for aer_idx = 0, n_elements(aerosol_element_names)-1 do begin
               if stregex(names[name_idx], aerosol_element_names[aer_idx], /BOOLEAN) then begin
                   if names[name_idx] NE curr_aer_name then begin
                       curr_aer_idx += 1
                       curr_aer_name = names[name_idx]
                   endif

                   names[name_idx] = string(curr_aer_idx, FORMAT=aerosol_group_format)
                   break
               endif
           endfor

       endif

       ;; Remove any extra _XXX in names
       while stregex(names[name_idx], '_[0-9]{3}$', /BOOLEAN) do begin
           names[name_idx] = strmid(names[name_idx], 0, strlen(names[name_idx])-4)
       endwhile

   endfor

   return, names

END
