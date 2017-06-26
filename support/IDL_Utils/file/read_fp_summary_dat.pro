function read_fp_summary_dat, file

   matrix_obj = Obj_New('MATRIX_FILE', file)

   ncases = matrix_obj->Get_Num_Rows()

   ;; get names
   tags = matrix_obj->Get_All_Column_Labels()

   command = 'elt_ = {'
   for t = 0, n_elements(tags)-1 do begin
       nam = strupcase(strtrim(tags[t]))
       this = '0.0D0'
       if strpos(nam, 'CASE') ne -1 then begin
           this = '" "'
           nam = 'CASE_NUMBER'  ; "case" is a reserved word in IDL
       endif
       if strpos(nam, 'YEAR') ne -1 then this = '0'
       if strpos(nam, 'DAY') ne -1 then this = '0'
       if strpos(nam, 'FRAME_TIME') ne -1 then this = '" "'
       if strpos(nam, 'ITERATIONS') ne -1 then this = '0'
       if strpos(nam, 'OUTCOME') ne -1 then this = '0'
       if strpos(nam, 'FLAG') ne -1 then this = '" "'
       
       command = command + nam + ':' + this
       if (t ne (n_elements(tags)-1)) then command = command + ', '
   endfor

   command = command + '}'
   err = execute(command)

   all = replicate(elt_, ncases)

   file_data = matrix_obj->Get_Data(/NOCASTDOUBLE)

   for case_idx = 0, ncases-1 do begin
       for col_idx = 0, n_elements(tags)-1 do begin
           if (file_data[col_idx, case_idx] EQ "N/A") then begin
               all[case_idx].(col_idx) = -999.0
           endif else begin
               all[case_idx].(col_idx) = file_data[col_idx, case_idx]
           endelse
       endfor
   endfor

   return, all

END
