function strswitch, change_str, from_str, to_str
   
   pos = strpos(change_str, from_str)
   while pos gt 0 do begin
       change_str = strmid(change_str, 0, pos) + to_str + strmid(change_str, pos + strlen(from_str))
       pos = strpos(change_str, from_str)
   endwhile

   return, change_str
end

function make_osim_log_struct, nameline, indexline

   x = strsplit(indexline, count=n, len=lens)
   b=0
   com = 'elt_ = {'
   for i = 0, n-2 do begin
       a = b+1
       b = x[i+1]+lens[i+1]-1
       nam = strtrim(strmid(nameline, a, b-a+1),1)
       ;; remove last parantheses if necessary
       p = strpos(nam, ')')
       if p ne -1 then nam = strtrim(strmid(nam, 0, p),0)  
       nam = strswitch(nam, ' (', '_')
       nam = strswitch(nam, '(', '_')
       nam = strswitch(nam, '^','') ; remove hats (in units)
       nam = strswitch(nam, '/','_') ; exchange / for _
       nam = strswitch(nam,' ','_') ; exchange blanks for underscores
       nam = strswitch(nam,'?', '') ; remove question marks
       
       this = '0.'
       lnam = strlowcase(nam)
       if strpos(lnam, 'frame') ne -1 then this = '0'
       if strpos(lnam, 'index') ne -1 then this ='0'
       if strpos(nam, 'ID') ne -1 then this = '0ULL'
       if strpos(lnam, 'sounding') ne -1 then this = '0'
       if strpos(lnam, 'type') ne -1 then this = '0'
       if strpos(lnam, 'modis') ne -1 then this = '0'
       if strpos(lnam, 'igbp') ne -1 then this = '0'
       com = com + nam + ':' + this
       if i ne (n-2) then com = com + ', '
   endfor
   com = com + '}'
   err = execute(com)

   return, elt_

end

function read_orbit_simulator_log, file

   if ~file_test(file) then begin
       print, 'Orbit simulator log file ' + file + ' not found.'
       RETURN, -1
   endif

   line = ''
   openr, lun, file, /get_lun
   
   repeat begin
       readf, lun, line
   endrep until strlen(line) GT 100
   nameline = line
   indexline = ''
   readf, lun, indexline
   readf, lun, line             ; last blank line before data begin
   
   elt_ = make_osim_log_struct(nameline, indexline)
   
   out = replicate(elt_, 20000) ; maximum number of profiles in a log file
   
   ;; now read the actual lines
   n = 0L
   while ~EOF(lun) do begin
       readf, lun, elt_
       out[n] = elt_
       n += 1
   endwhile
   close, lun
   free_lun, lun
   out = out[0:n-1]
   
   return, out

END
