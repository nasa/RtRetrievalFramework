pro Plot_Spectrum_Residuals, case_name, run_dir, plot_dir, overwrite

   errStat = 0
   ;Catch, errStat
   IF (errStat NE 0) THEN BEGIN
       ON_ERROR, 2
       Catch, /CANCEL
       Message, /CONTINUE, !Error_State.msg
       Message, /CONTINUE, 'Error encountered can not continue'
       
       if n_elements(device_init) gt 0 then $
         device, /close_file

       if n_elements(charsize_orig) gt 0 then $
         !p.charsize = charsize_orig

       if n_elements(multi_orig) gt 0 then $
         !p.multi = multi_orig
       
       if n_elements(x_style_orig) gt 0 then $
         !x.style = x_style_orig

       if n_elements(y_style_orig) gt 0 then $
         !y.style = y_style_orig
       
       return
   ENDIF

   if n_elements(plot_dir) eq 0 then $
     plot_dir = './'

   if n_elements(overwrite) eq 0 then $
     overwrite = 1

   ps_filename = plot_dir + '/' + case_name + '_spec_residuals.ps'

   existing = file_search(ps_filename)
   if existing[0] ne '' and not overwrite then begin
       print, 'Skipping existing: ', ps_filename
       return
   endif

   ;; Set up for plot parts
   leg_names = ['Measured', 'Calculated']
   txt_colors = [1, 1]
   item_colors = [1, 2]
   item_psym = [0, 0]
   ;item_lines = [0, 2]
   item_lines = [0, 0]

   charsize_orig = !p.charsize
   !p.charsize = 0.5

   legend_spacing = 0
               
   ;; Search for per iteration file
   fnames_conv  = file_search(run_dir + '/out/rad_conv.dat.iter*')
   fnames_meas  = file_search(run_dir + '/out/rad_meas.dat.iter*')
   n_files = size(fnames_conv, /dimension)

   ;; If no per iteration files, use last iteration
   if n_files[0] eq 0 then begin
       fnames_conv  = file_search(run_dir + '/out/rad_conv.dat')
       fnames_meas  = file_search(run_dir + '/out/rad_meas.dat')
       n_files = size(fnames_conv, /dimension)
   endif

   if n_files[0] eq 0 then $
     Message, 'Could not find any rad_conv.dat files for run dir: ' + run_dir

   print, 'Creating: ', ps_filename

   set_plot, 'ps'
   device, /color, bits_per_pixel=8, filename=ps_filename, /portrait
   device_init = 1
   for numfile = 0, n_files[0]-1 do begin
       fname_out  = fnames_conv[numfile]
       fname_in   = fnames_meas[numfile]

       plot_title = case_name + ' Spectrum Residuals'
       iter_pos = strpos(fname_out, 'iter')
       if iter_pos ge 0 then $
         plot_title = plot_title + ' Iteration = ' + strmid(fname_out, iter_pos + 4) + ', '

       ;; read simulated measurement spectrum
       rad_obj = obj_new('matrix_file', fname_in)

       start_pixels_str = rad_obj->Get_Header_Keyword('Start_Pixels')
       pix_parts = strsplit(start_pixels_str, /extract)
       num_windows = n_elements(pix_parts)

       Num_Cols = rad_obj->get_num_columns()
       Num_Rows = rad_obj->get_num_rows()

       if Num_Cols eq 0 or Num_Rows eq 0 then $
         Message, 'Measured radiance file is empty: ' + fname_in

       in_rad = rad_obj->get_data()
       obj_destroy, rad_obj

       ;; read retrieved spectrum
       rad_obj = obj_new('matrix_file', fname_out)
       out_rad = rad_obj->get_data()
       obj_destroy, rad_obj

       ;; compute residual spectrum
       residual = out_rad[*,*]
       residual[2,*] = in_rad[2,*]-out_rad[2,*]

       for win_index = 0, num_windows - 1 do begin
           plot_title_win = plot_title + ' Window = ' + string(win_index, format='(I0)')

           multi_orig = !p.multi
           x_style_orig = !x.style
           y_style_orig = !y.style

           !p.multi = [0, 0, 2, 0, 0]
           !x.style = 1
           !y.style = 1

           Start_1 = pix_parts[win_index]
           if win_index lt num_windows-1 then begin
               Start_2 = pix_parts[win_index+1]
           endif else begin
               Start_2 = Num_Rows + 1
           endelse

           wavelength = in_rad[0, Start_1 - 1 : Start_2 - 2]
           wavenumber = in_rad[1, Start_1 - 1 : Start_2 - 2]
           in_spec = in_rad[2, Start_1 - 1 : Start_2 - 2]
           out_spec = out_rad[2, Start_1 - 1 : Start_2 - 2]
           resd_spec = 100. * (residual[2, Start_1 - 1 : Start_2 - 2] / max(in_spec))
           smax = max([max(in_spec), max(out_spec)])
           rmax = max([abs(min(resd_spec)), max(resd_spec)])

           plot, wavelength, in_spec, /nodata, color=1, position=[0.13, 0.10, 0.95, 0.82], $
;                 xrange = [min(wavenumber), max(wavenumber)], xtitle = 'Wavenumber (cm!u-1!n)', $
                 xrange = [min(wavelength), max(wavelength)], xtitle = 'Wavelength (micron)', $
                 xtickname =['', '', '', '', '', '', '', '', '', ''], yrange = [0., smax], $
                 ytitle = 'Radiance'     
           oplot, wavelength, in_spec, linestyle=item_lines[0], psym=item_psym[0], color = item_colors[0]
           oplot, wavelength, out_spec, linestyle=item_lines[1], psym=item_psym[1], color = item_colors[1]

           Legend_DR, leg_names, TEXTCOLORS=txt_colors, COLORS=item_colors, PSYM=item_psym, $
                      LINESTYLE=item_lines, CORNERS=corners, POS=[2,2], SPACING=legend_spacing

           xydims = [corners[2]-corners[0],corners[3]-corners[1]]
           chdim=[!d.x_ch_size/float(!d.x_size),!d.y_ch_size/float(!d.y_size)]
           pos = [!x.window[1]-chdim[0]/2.0-xydims[0], !y.window[0]+xydims[1]/2.0+chdim[1]]

           polyx  = [ pos[0], pos[0] + xydims[0], pos[0] + xydims[0], pos[0] ]
           polyy  = [ pos[1], pos[1], pos[1] - xydims[1], pos[1] - xydims[1] ]
           
           ;; Fill the area where the legend will be drawn so that no data will
           ;; leak through and obscure the legend
           polyfill, polyx, polyy, /NORM, COLOR=0

           Legend_DR, leg_names, TEXTCOLORS=txt_colors, COLORS=item_colors, PSYM=item_psym, $
                      LINESTYLE=item_lines, POS=pos, SPACING=legend_spacing


           plot, wavelength, resd_spec, /nodata, color = 1, title=plot_title_win, $
                 position = [0.13, 0.83, 0.95, 0.95], xtitle = ' ', $
                 xtickname =[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], $
                 ytitle = 'Residual (%)', yrange = [-1*rmax, rmax]

           oplot, wavelength, dblarr(n_elements(resd_spec)), linestyle = 0, color = 1 ; zero line
           oplot, wavelength, resd_spec, linestyle = 0, color = 1

           !p.multi = multi_orig
           !x.style = x_style_orig
           !y.style = y_style_orig
       endfor

   endfor

   !p.charsize = charsize_orig

   device, /close_file

end

