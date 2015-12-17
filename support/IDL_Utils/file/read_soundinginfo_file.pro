FUNCTION Read_Soundinginfo_File, soundinginfo_file, DATA_NAMES=data_names

   ;; Number of values read from the file
   num_vals_to_read = 40

   ;; Read soundinginfo file
   OpenR, lun, soundinginfo_file, /GET_LUN

   data_names  = StrArr( num_vals_to_read )
   data_values = StrArr( num_vals_to_read )

   val_count = 0
   line_str = ''
   while (eof(lun) ne 1 and val_count lt num_vals_to_read) do begin
       Readf, lun, line_str

       if strpos(strtrim(line_str,2), '#') ne 0 and strpos(line_str, '=') ge 0 then begin
           line_parts = strsplit(line_str, '=', /EXTRACT)
           
           curr_name = StrTrim(line_parts[0], 2)      
           curr_val  = StrTrim(line_parts[1], 2)
           
           data_names [ val_count ] = curr_name
           data_values[ val_count ] = curr_val

           val_count = val_count + 1
       endif
   endwhile

   data_names = data_names [0:val_count-1]
   data_vals  = data_values[0:val_count-1]

   Free_Lun, lun
   
   return, data_vals

END
