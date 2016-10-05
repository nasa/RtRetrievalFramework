PRO Draw_Info_Box, text_lines, POS=pos, CORNERS=corners, OUTLINE=do_outline, _EXTRA=extra

   ;; Get number of columns and rows to help with sizing the bounding box
   num_char_cols = Max(strlen(text_lines))
   num_char_rows = N_Elements(text_lines)
   
   chdim = [!d.x_ch_size/float(!d.x_size),!d.y_ch_size/float(!d.y_size)]

   ;; Use default poistion if not specified
   IF N_Elements(pos) EQ 0 THEN $
     pos = [ !x.window[1]-(chdim[0]/2.0)*(num_char_cols+1), !y.window[1]-(chdim[1]/2.0)* double(num_char_rows+2) ] ;; Top right
   
   ;; Figure out the corners of the box to bound the region
   polyx  = [ pos[0]+chdim[0]/4.0, pos[0]+(chdim[0]/2.0)*num_char_cols, pos[0]+(chdim[0]/2.0)*num_char_cols, pos[0]+chdim[0]/4.0 ]
   polyy  = [ pos[1]+(chdim[1]/2.0)*double(num_char_rows+1), pos[1]+(chdim[1]/2.0)*double(num_char_rows+1), pos[1]+(chdim[1]/2.5), pos[1]+(chdim[1]/2.5) ]

   ;; Return this array to indicate where the corners are
   corners = [ polyx[0], polyy[0], polyx[1], polyy[1] ]
   
   ;; Fill the area where the box will be drawn so that no data will
   ;; leak through and obscure the information
   polyfill, polyx, polyy, /NORM, COLOR=0
   
   ;; Draw text lines
   for line_index = 0, num_char_rows-1 do begin
       y_mult = num_char_rows - line_index

       xyouts, pos[0]+chdim[0]/2.0, pos[1]+(chdim[1]/2.0) * double(y_mult), $
               text_lines[line_index], /NORMAL, _EXTRA=extra
   endfor

   ;; Draw outline of box
   IF Keyword_Set(do_outline) THEN BEGIN
       boxx = replicate(0.0D, N_Elements(polyx) + 1)
       boxy = replicate(0.0D, N_Elements(polyy) + 1)
   
       boxx[0:N_Elements(polyx)-1] = polyx
       boxx[N_Elements(polyx)] = polyx[0]
       
       boxy[0:N_Elements(polyy)-1] = polyy
       boxy[N_Elements(polyy)] = polyy[0]
   
       PlotS, boxx, boxy, /NORM, COLOR=1
   ENDIF

END
