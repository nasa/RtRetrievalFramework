FUNCTION Load_True_Data, run_dir, true_type, HAS_TRUE=has_true, pressure_profile=pressures

   has_true = 0B

   true_base_dir = run_dir + '/true'
   true_files_map = [ [ 'CO2',        'atmosphere*.dat', 'CO2'     ,    'APPEND' ], $
                      [ 'H2O',        'atmosphere*.dat', 'H2O'     ,    'APPEND' ], $
                      [ 'TEMP',       'atmosphere*.dat', 'T'       ,    'APPEND' ], $
                      [ 'PSURF',      'psurf*.dat',      'PSURF'   ,    'APPEND' ], $
                      [ 'AEROSOL',    'aerosol*.dat',    '.*'      ,    'ADD'    ], $
                      [ 'CONT',       'aerosol*.dat',    'CONT'    ,    'ADD' ], $
                      [ 'OCEANIC',    'aerosol*.dat',    'OCEANIC' ,    'ADD' ], $
                      [ 'AERO1',      'aerosol*.dat',    'CLAT_0'  ,    'ADD' ], $
                      [ 'AERO2',      'aerosol*.dat',    'CLAT_1'  ,    'ADD' ], $
                      [ 'AERO1',      'aerosol*.dat',    'AERO1'   ,    'ADD' ], $
                      [ 'AERO2',      'aerosol*.dat',    'AERO2'   ,    'ADD' ], $
                      [ 'WC',         'aerosol*.dat',    'WC_?.*'  ,    'ADD' ], $
                      [ 'IC',         'aerosol*.dat',    'IC_?.*'  ,    'ADD' ], $
                      [ 'LOW',        'aerosol*.dat',    'LOW.?_'  ,    'ADD' ], $
                      [ 'MID',        'aerosol*.dat',    'MID.?_'  ,    'ADD' ], $
                      [ 'HIGH',       'aerosol*.dat',    'HIGH.?_' ,    'ADD' ], $
                      [ 'ALBEDO',     'albedo*.dat',     'ALBEDO_.',    'APPEND' ], $
                      [ 'WINDSPEED',  'windspeed*.dat',  'WINDSPEED_.', 'APPEND' ] $
                    ]

   ; Where source type data is loaded from
   map_where = where( true_type eq true_files_map[0, *] )
   
   if map_where[0] EQ -1 then begin
       Message, /CONTINUE, 'Type ' + true_type + ' not found in true file map'
       Return, -1
   endif

   has_true = 0B

   for map_idx = 0, n_elements(map_where)-1 do begin
       true_data_file = file_search(true_base_dir + '/' + true_files_map[1, map_where[map_idx]]) 
       if strlen(true_data_file[0]) eq 0 then begin
           Message, /CONTINUE, 'True data file not found as: ' + true_base_dir + '/' + true_files_map[1, map_where[map_idx]]
           Return, -1
       endif
      
       file_column_name = true_files_map[2, map_where[map_idx]]
       add_columns      = true_files_map[3, map_where[map_idx]] eq 'ADD'

       true_data_obj = Obj_New('MATRIX_FILE', true_data_file[0])
       data_columns = true_data_obj->get_all_column_labels()

       for col_idx = 0, N_Elements(data_columns)-1 do begin
           col_data = true_data_obj->get_column_data(col_idx)

           if stregex(data_columns[col_idx], 'Pressure', /BOOLEAN, /FOLD_CASE) then begin
               pressures = col_data[*]
               continue
           endif

           col_match = stregex(data_columns[col_idx], file_column_name, /BOOLEAN, /FOLD_CASE)

           if col_match then begin
               if n_elements(true_data) gt 0 then begin
                   if add_columns then begin
                       true_data = true_data + col_data[*]
                   endif else begin
                       new_true = DblArr( N_Elements(true_data) + N_Elements(col_data) )
                       new_true[0:N_Elements(true_data)-1] = true_data
                       new_true[N_Elements(true_data):*] = col_data[*]
                       true_data = new_true
                   endelse
               endif else begin
                   true_data = col_data[*]
               endelse

               has_true = 1B
           endif
       endfor
   endfor

   if n_elements(true_data) eq 0 then $
     return, -1 $
   else $
     return, true_data
END

FUNCTION Get_Type_Size, type_name, sv_names
   where_type = where(sv_names eq type_name, type_size)
   return, type_size
END

FUNCTION Get_Type_SV_Range, type_name, sv_len, sv_names
   where_type = where(sv_names eq type_name)

   return, [where_type[0], where_type[n_elements(where_type)-1]]
END

FUNCTION Fix_SV_Names_Size, sv_names, desired_size
   LEVEL_ITEMS = [ 'CO2', $
                   'H2O', $
                   'TEMP',$
                   'AEROSOL', $
                   'CONT',$
                   'OCEANIC',$
                   'AERO1',$
                   'AERO2',$
                   'WC',$
                   'IC'$
                 ]

   MAX_ATTEMPT = 20

   if n_elements(sv_names) eq desired_size then $
     return, sv_names

   ; How do we add leveled values to try
   if n_elements(sv_names) lt desired_size then begin
       direction = 1 
   endif else begin
       direction = -1
   endelse

   sv_types = sv_names[ uniq(sv_names) ]
   
   current_size = n_elements(sv_names)
   attempt = 0
   new_names = strarr(desired_size)
   while (attempt lt MAX_ATTEMPT) do begin
       if direction eq 1 and current_size ge desired_size then $
         break $
       else if direction eq -1 and current_size le desired_size then $
         break

       new_names[*] = ''
       attempt = attempt+1

       new_sv_idx = 0
       for curr_type = 0, n_elements(sv_types)-1 do begin
           where_type_in_sv = where(sv_names eq sv_types[curr_type], type_size)
           where_type_in_lvl = where(LEVEL_ITEMS eq sv_types[curr_type], lv_size)
           
           copy_size = 0
           if lv_size gt 0 and type_size gt 1 then begin
               copy_size = type_size + (attempt*direction)
           endif else begin
               copy_size = type_size
           endelse

           new_names[new_sv_idx:new_sv_idx+copy_size-1] = sv_types[curr_type]

           new_sv_idx = new_sv_idx+copy_size
       endfor
       
       current_size = new_sv_idx
   endwhile

   if (current_size ne desired_size) then begin
       print, 'Could not adjust statevector names to desired size:', desired_size
       return, ['']
   endif   

   return, new_names[0:desired_size-1]
END

PRO Plot_StateVector_Profiles, case_name, run_dir, plot_dir, overwrite

   !except = 2

   if n_elements(plot_dir) eq 0 then $
     plot_dir = './'

   if n_elements(overwrite) eq 0 then $
     overwrite = 1

   non_pres_types = [ 'ALBEDO', 'DISP' ]

   svn_rev_str = ''
   ;svn_rev_str = 'SVN Revision: ' + StrMid(run_dir, StrPos(run_dir, '-', /reverse_search) + 1 )

   ;;;;
   ;; File control paramters
   plot_file_format = '("' + plot_dir + '/' + case_name +'_",A0,"_sv.ps")'

   pres_filename = file_search(run_dir + '/in/scene/atmosphere/pressure*.dat')

   ;; Check for an atmosphere file
   if pres_filename[0] eq "" then $
     pres_filename = file_search(run_dir + '/in/scene/atmosphere/atmosphere*.dat')

   ;; Check for an output atmosphere file
   if pres_filename[0] eq "" then $
     pres_filename = file_search(run_dir + '/out/atmosphere.dat*.iter00')

   if pres_filename[0] eq "" then $
     pres_filename = file_search(run_dir + '/out/atmosphere.dat.spec01.iter00')

   if pres_filename[0] eq "" then begin
       Message,'Pressure file could not be found for run directory: ' + run_dir, /CONTINUE
       Return
   endif

   ;; Try first for per iteration state vector files
   sv_filenames = file_search(run_dir + '/out/statevector.dat.iter??')
   sv_apriori_col  = 1
   sv_solution_col = 2
   sv_update_col   = 3

   ;; Try now just for final state vector file
   if sv_filenames[0] eq "" then begin
       sv_filenames = file_search(run_dir + '/out/statevector.dat')
       sv_apriori_col  = 1
       sv_solution_col = 2
       sv_update_col   = 3
   endif

   ;; Use error analysis state vector file
   if sv_filenames[0] eq "" then begin
       sv_filenames = file_search(run_dir + '/out/control1/statevector.dat')
       sv_apriori_col  = 0
       sv_solution_col = 1
   endif

   ;; Still could not find state vector so die!
   if sv_filenames[0] eq "" then begin
       Message, 'State vector file cound not be found for run directory: ' + run_dir, /CONTINUE
       Return
   endif

   sv_names_file = file_search(run_dir + '/out/aggregator/sv_names.dat')
   if sv_names_file[0] eq "" then begin
       Message, 'Statevector names file cound not be found for run directory: ' + run_dir, /CONTINUE
       Return
   endif

 
   ;; Plotting parameters
   make_ps = 1

   if n_elements(sv_filenames) gt 1 then begin
       num_plot_cols = 2
       num_plot_rows = 2
   endif else begin
       num_plot_cols = 1
       num_plot_rows = 1
   endelse

   plot_multi_array = [0, num_plot_cols, num_plot_rows, 0, 0]
   plots_per_page = num_plot_cols * num_plot_rows

   charsize_orig = !p.charsize

   !p.charsize = 0.5
   legend_spacing = 0

   item_lines  = [0, 2, 7]
   item_colors = [1, 2, 7]

   title_format='("'+case_name+' ",A0," State Vector Comparison")'

   pressure_y_title_str='Pressure (hPa)'
   generic_y_title_str='Index'

   ;;;;;;;;;;;;;;;;;

   sv_names_all = Read_Statevector_Names(sv_names_file)
   sv_names_aer = Read_Statevector_Names(sv_names_file, /group_aerosol)

   ;; Calculate  unique state vector element types
   sv_types_all = sv_names_all[ uniq(sv_names_all) ]
   sv_types_aer = sv_names_aer[ uniq(sv_names_aer) ]
  
   where_aero_type = where( stregex(sv_types_aer, 'AEROSOL_', /BOOLEAN, /FOLD_CASE) eq 1 )
   ;; Push aerosol down into one type
   if where_aero_type[0] ne -1 then begin
       aerosol_types = sv_types_all[where_aero_type]
       sv_types_all  = [sv_types_all, 'AEROSOL']
   endif

   ;; Load true PSURF if it exits
   true_search_result = Load_True_Data( run_dir, 'PSURF', has_true=has_true_psurf )
   if has_true_psurf then $
     true_psurf = true_search_result / 100.0D

   ;; Load pressure data
   pres_obj = Obj_New('MATRIX_FILE', pres_filename[0])
   pres_data_full = pres_obj->get_data() / 100.0D
   if (size(pres_data_full))[0] gt 1 then $
     pres_data_full = pres_data_full[0,*]

   pressure_range = [pres_data_full[n_elements(pres_data_full)-1], pres_data_full[0]]

   if n_elements(true_psurf) gt 0 then begin
       pressure_range[0] = true_psurf
   endif

   ;; Loop over all state vector types
   for type_idx = 0, N_Elements(sv_types_all)-1 do begin
       curr_sv_type = sv_types_all[type_idx]

       ;; Skip PSURF element, instead write text in box
       if curr_sv_type eq 'PSURF' then continue

       ;; Find the location in the names where the type exists
       if curr_sv_type eq 'AEROSOL' then begin
           type_srch_name = aerosol_types[0]
       endif else begin
           type_srch_name = curr_sv_type
       endelse

       ;; Set up plot file for this type
       if make_ps then begin
           ps_filename=String(curr_sv_type,FORMAT=plot_file_format)

           existing = file_search(ps_filename)
           if existing[0] ne '' and not overwrite then begin
               print, 'Skipping existing: ', ps_filename
               continue
           endif

           entry_device = !d.name
           set_plot, 'PS'

           orig_multi = !p.multi
           !p.multi = plot_multi_array

           print, 'Creating: ', ps_filename
           device, /color, bits_per_pixel=8, filename=ps_filename, /portrait, /helvetica
       endif

       ;; Load true data for type
       true_pressures = -1
       has_sv_true = 0B
       true_search_result = Load_True_Data( run_dir, curr_sv_type, HAS_TRUE=has_sv_true, pressure_profile=true_pressures )

       if has_sv_true then $
         sv_true = true_search_result

       ;; Initialize min/max from true data if it exists
       if has_sv_true then begin
           data_min = min(sv_true)
           data_max = max(sv_true)
       endif else begin
           data_min = (machar(/double)).xmax
           data_max = (machar(/double)).xmin
       endelse

       ;; Load state vector file data
       max_type_size = 100L ;Get_Type_Size( type_srch_name, sv_names_curr )
       sv_solution_matrix = replicate(0D, n_elements(sv_filenames), max_type_size)
       sv_apriori_matrix = replicate(0D, n_elements(sv_filenames), max_type_size)
       sv_error_matrix = replicate(0D, n_elements(sv_filenames), max_type_size)
       sv_type_iter_sizes = replicate(0D, n_elements(sv_filenames))

       ;; Find location of PSURF for adjusting pressure
       check_psurf = Where(sv_types_all eq 'PSURF')
       if check_psurf[0] ne -1 then begin
           sv_psurfs = replicate(0D, n_elements(sv_filenames))
       endif

       ;; False that aer in log until otherwise known
       aer_in_log = 0B

       labels = file_basename(sv_filenames)
       for l_idx = 0, n_elements(labels)-1 do begin
           iter_pos = strpos(labels[l_idx], 'iter')
           if iter_pos ge 0 then $
             labels[l_idx] = 'iteration ' + strmid(labels[l_idx], iter_pos + 4)$
           else $
             labels[l_idx] = 'solution'
       endfor

       ;; Loop over sv filenames and aggregate data and get max/mins
       for svFileIdx = 0, N_Elements(sv_filenames)-1 do begin
           sv_obj = Obj_New('MATRIX_FILE', sv_filenames[svFileIdx])
           sv_len = sv_obj->get_num_rows()

           sv_names_curr = fix_sv_names_size(sv_names_all, sv_len)

           if sv_names_curr[0] eq '' then begin
               return
           endif

           type_range = Get_Type_SV_Range(type_srch_name, sv_len, sv_names_curr)

           sv_beg = type_range[0]
           sv_end = type_range[1]

           sv_type_iter_sizes[svFileIdx] = sv_end-sv_beg+1

           sv_apriori    = (sv_obj->get_column_data( sv_apriori_col ))[sv_beg:sv_end]
           sv_solution   = (sv_obj->get_column_data( sv_solution_col ))[sv_beg:sv_end]

           if n_elements(sv_update_col) gt 0 then begin
               sv_solution   = sv_solution + (sv_obj->get_column_data( sv_update_col ))[sv_beg:sv_end]
           endif

           ;; Add all other aerosol types to the current one
           where_in_aer_types = where(type_srch_name eq aerosol_types)
           if where_in_aer_types[0] ne -1 then begin
               where_neg = where(sv_solution lt 0.0, neg_count)
               aer_in_log = neg_count gt 0
               if aer_in_log then begin
                   sv_apriori  = exp(sv_apriori)
                   sv_solution = exp(sv_solution)
               endif

               if has_sv_true then begin
                   where_neg = where(sv_true lt 0.0, neg_count)
                   true_in_log = neg_count gt 0
                   if true_in_log then $
                     sv_true = exp(sv_true)
               endif
               
               if sv_types_all[type_idx] eq 'AEROSOL' and n_elements(aerosol_types) gt 1 then begin

                   for aer_idx = 1, n_elements(aerosol_types)-1 do begin
                       
                       aer_apriori_add  = (sv_obj->get_column_data( sv_apriori_col ))[sv_beg:sv_end]
                       aer_solution_add = (sv_obj->get_column_data( sv_solution_col ))[sv_beg:sv_end]

                       if aer_in_log then begin
                           aer_apriori_add  = exp(aer_apriori_add)
                           aer_solution_add = exp(aer_solution_add)
                       endif

                       sv_apriori  = sv_apriori  + aer_apriori_add
                       sv_solution = sv_solution + aer_solution_add
                   endfor
               endif
           endif
       
           if sv_obj->Get_Num_Columns() ge 5 then begin
               sv_error = (sv_obj->get_column_data( 4 ))[sv_beg:sv_end]

               if aer_in_log then $
                   sv_error = exp(sv_error)

               sv_error_matrix(svFileIdx, 0:sv_type_iter_sizes[svFileIdx]-1)  = sv_error
           endif

           if check_psurf[0] ne -1 then begin
               psurf_range = Get_Type_SV_Range('PSURF', sv_len, sv_names_curr)               

               sv_psurfs[svFileIdx] = (sv_obj->get_column_data( sv_solution_col ))[psurf_range[0]] / 100.0D
               pressure_range[0] = max([sv_psurfs[svFileIdx], pressure_range[0]])
           endif

           sv_solution_matrix(svFileIdx, 0:sv_type_iter_sizes[svFileIdx]-1) = sv_solution
           sv_apriori_matrix(svFileIdx, 0:sv_type_iter_sizes[svFileIdx]-1)  = sv_apriori
       endfor

       title_str = String(curr_sv_type, FORMAT=title_format)
       x_title_str = curr_sv_type

       ;; Loop over iterations again and create plots
       for iter_idx = 0, N_Elements(labels)-1 do begin
           sv_solution = sv_solution_matrix[iter_idx, 0:sv_type_iter_sizes[iter_idx]-1]
           sv_apriori  = sv_apriori_matrix[iter_idx, 0:sv_type_iter_sizes[iter_idx]-1]
           sv_error    = sv_error_matrix[iter_idx, 0:sv_type_iter_sizes[iter_idx]-1]
           pres_data   = pres_data_full[0:sv_type_iter_sizes[iter_idx]-1]

           if n_elements(sv_psurfs) gt 0 then begin
               ret_psurf = sv_psurfs[iter_idx]
               pres_data[N_Elements(pres_data)-1] = ret_psurf
           endif

           where_non_pres_type = where(non_pres_types eq curr_sv_type)
           if where_non_pres_type eq -1 and sv_type_iter_sizes[iter_idx] gt 1 then begin
               y_title_str = pressure_y_title_str
               y_axis = pres_data
               y_range_array = pressure_range

               item_psym   = [-1, -5, -6]
               p_symsize = 0.5
           endif else begin
               y_title_str = generic_y_title_str

               ;if (sv_type_iter_sizes[iter_idx]) gt 1 then $
               y_axis = indgen(sv_type_iter_sizes[iter_idx])
               ;else $
               ;  y_axis = indgen(3) - 1
               
               y_range_array = [y_axis[N_Elements(y_axis)-1]+1, y_axis[0]-1]

               item_psym   = [1, 5, 6]
               p_symsize = 0.8
           endelse
           
           data_min = min([min(sv_solution), min(sv_apriori)])
           data_max = max([max(sv_solution), max(sv_apriori)])

           if has_sv_true then begin
               data_min = min([data_min, min(sv_true)])
               data_max = max([data_max, max(sv_true)])
           endif
           
           data_min = data_min - (data_max-data_min)*.1
           data_max = data_max + (data_max-data_min)*.1

           plot, sv_apriori, y_axis, linestyle=item_lines[0], color=item_colors[0], psym=item_psym[0], title=title_str, xtitle=x_title_str, ytitle=y_title_str, xrange=[data_min, data_max], yrange=y_range_array, symsize=p_symsize, xstyle=1

           if n_elements(sv_error) gt 0 then begin
               if aer_in_log then begin

                   err_plot, y_axis, $
                             sv_apriori + (sv_apriori * sv_error) / 2.0, $
                             sv_apriori - (sv_apriori * sv_error) / 2.0, $
                             /flip, xrange=[data_min, data_max], yrange=[y_axis[N_Elements(y_axis)-1], y_axis[0]], color=1
               endif else begin
                   err_plot, y_axis, $
                             sv_apriori - sv_error / 2.0, $
                             sv_apriori + sv_error / 2.0, $
                             /flip, xrange=[data_min, data_max], yrange=[y_axis[N_Elements(y_axis)-1], y_axis[0]], color=1
               endelse
           endif

           if has_sv_true then begin
               if true_pressures[0] eq -1 then $
                 true_axis = y_axis $
               else $
                 true_axis = true_pressures / 100.0D0

               oplot, sv_true, true_axis, linestyle=item_lines[1], color=item_colors[1], psym=item_psym[1], symsize=p_symsize
           endif

           oplot, sv_solution, y_axis, linestyle=item_lines[2], color=item_colors[2], psym=item_psym[2], symsize=p_symsize
           
           legNames  = ['apriori', 'true', labels[iter_idx]]
           txtColors = [1, 1, 1]
           
           Legend_DR, legNames, TEXTCOLORS=txtColors, COLORS=item_colors, PSYM=item_psym, LINESTYLE=item_lines, CORNERS=corners, POS=[2,2], SPACING=legend_spacing

           xydims = [corners[2]-corners[0],corners[3]-corners[1]]
           chdim=[!d.x_ch_size/float(!d.x_size),!d.y_ch_size/float(!d.y_size)]

           ; Top middle-ish
           pos = [!x.window[0]+(!x.window[1]-!x.window[0])/2.0+(xydims[1]-xydims[0])/2.0, $
                  !y.window[1]-(chdim[1]/2.0) ]

           ; Bottom right
           ;pos = [!x.window[1]-chdim[0]/2.0-xydims[0], !y.window[0]+xydims[1]/2.0+chdim[1]]

           polyx  = [ pos[0], pos[0] + xydims[0], pos[0] + xydims[0], pos[0] ]
           polyy  = [ pos[1], pos[1], pos[1] - xydims[1], pos[1] - xydims[1] ]
           
           ;; Fill the area where the legend will be drawn so that no data will
           ;; leak through and obscure the legend
           ;polyfill, polyx, polyy, /NORM, COLOR=0

           Legend_DR, legNames, TEXTCOLORS=txtColors, COLORS=item_colors, PSYM=item_psym, LINESTYLE=item_lines, POS=pos, SPACING=legend_spacing

           if (strlen(svn_rev_str) gt 0) and (iter_idx mod plots_per_page) eq 0 then $
             xyouts, 1.0-chdim[0]/2.0, 0.0+chdim[1]/2.0, svn_rev_str, alignment=1.0, /normal, color=1


           ;; Draw info box lines
           info_box_lines = "-1"
           
           if N_Elements(ret_psurf) gt 0 then begin

               psurf_info_lines = strarr(3)

               psurf_info_lines[0] = 'PSURF:'
               psurf_info_lines[1] = string(ret_psurf, '("Retr: ",F0.2," hPa")')

               if N_Elements(true_psurf) gt 0 then begin
                   psurf_info_lines[2] = string(true_psurf, '("True: ",F0.2," hPa")') 
               endif

               if info_box_lines[0] ne "-1" then begin
                   num_new_info_lines =  n_elements(info_box_lines) + n_elements(psurf_info_lines) 
                   new_info_lines = strarr(num_new_info_lines)
               
                   new_info_lines[0:n_elements(info_box_lines)-1] = info_box_lines
                   new_info_lines[n_elements(info_box_lines):*] = psurf_info_lines
                   info_box_lines = new_info_lines
               endif else begin
                   info_box_lines = psurf_info_lines
               endelse
                 
           endif

           if curr_sv_type eq 'CO2' then begin

               ret_pres = pres_data * 100.0D

               if N_Elements(ret_psurf) gt 0 then begin
                   x_target = xtarget_p_weight(ret_pres, sv_solution, psurf=ret_psurf*100.0D)
               endif else begin
                   x_target = xtarget_p_weight(ret_pres, sv_solution)
               endelse
               
               xco2_info_lines = strarr(3) ;; First line is blank to leave space in info box
               xco2_info_lines[1] = 'XCO2'
               xco2_info_lines[2] = string(x_target * 1e6, '("Retr: ",F0.3)')

               if info_box_lines[0] ne "-1" then begin
                   num_new_info_lines =  n_elements(info_box_lines) + n_elements(xco2_info_lines) 
                   new_info_lines = strarr(num_new_info_lines)
               
                   new_info_lines[0:n_elements(info_box_lines)-1] = info_box_lines
                   new_info_lines[n_elements(info_box_lines):*] = xco2_info_lines
                   info_box_lines = new_info_lines
               endif else begin
                   info_box_lines = xco2_info_lines
               endelse

           endif

           if info_box_lines[0] ne "-1" then $
             Draw_Info_Box, info_box_lines, /OUTLINE, COLOR=1

       endfor

       if make_ps then begin
           device, /close_file
           set_plot, entry_device

           !p.multi = orig_multi
       endif

   endfor

   !p.charsize = charsize_orig

END ;; pro plot_statevector_profiles
