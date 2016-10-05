FUNCTION Read_Out_Info_File, out_info_file, COLUMNS=col_names

   ;; Read out_info file
   OpenR, lun1, out_info_file, /GET_LUN

   while not eof(lun1) do begin
       line_str = ''
       ReadF, lun1, line_str

       if strpos(line_str, 'Additional Output file') ge 0 then begin
           file_id_str = line_str 
       endif else if strpos(line_str, 'kat kct') ge 0 then begin
           col_names_str = line_str
       endif else if strlen(line_str) gt 0 then begin
           data_str = line_str
       endif

   endwhile

   Free_Lun, lun1

   if n_elements(data_str) eq 0 then begin
       return, -1
   endif

   if data_str eq '' then begin
       return, -1
   endif

   ;; Fix data values so that last column has two error flags
   data_parts = StrSplit(data_str, /EXTRACT)

   data_values = StrArr(N_Elements(data_parts)-1)
   data_values[0:N_Elements(data_values)-2] = data_parts[0:N_Elements(data_values)-2]
   data_values[N_Elements(data_values)-1] = $
     StrJoin(data_parts[N_Elements(data_parts)-2:N_Elements(data_parts)-1], ', ')
   
   ;; Adjust column names taking into account number of spectrum
   col_name_parts = StrSplit(col_names_str, /EXTRACT)

   num_spec = (N_Elements(data_values) - 9) / 2

   ;; If this number does not make sense then might as well give up on
   ;; the file
   if num_spec le 0 then $
     return, -1
   
   col_names = StrArr(N_Elements(data_values))

   col_names[0:2] = col_name_parts[0:2]
   col_names[3:3+num_spec-1] = strmid(col_name_parts[2], 0, strpos(col_name_parts[2], '('))
   col_names[3+num_spec:3+(num_spec*2)-1] = strmid(col_name_parts[3], 0, strpos(col_name_parts[3], '('))
   col_names[3+(num_spec*2):*] =  col_name_parts[4:*]
               
   Return, data_values

END

