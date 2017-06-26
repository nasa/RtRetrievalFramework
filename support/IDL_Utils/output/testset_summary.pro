;####
;### DEPRECATED by Python routines
;###

function Format_Summary_Value, value
   double_format  = '(E0)'
   integer_format = '(I0)'

   if Valid_Num(value, /INTEGER) then begin
       Return, String( value, FORMAT=integer_format )
   endif else if Valid_Num(value) then begin
       Return, String( value, FORMAT=double_format )
   endif else begin
       Return, String( value )
   endelse
end

function Total_AOD_Values, layer_aod, indexes
   nlayers = n_elements(layer_aod)

   good_values = replicate(1, n_elements(indexes))

   outta_bound = where(indexes eq nlayers)
   if outta_bound ne -1 then $
     good_values[outta_bound] = 0

   where_good = where(good_values eq 1, num_good)

   if num_good gt 0 then begin
       indexes = indexes[where_good]
       return, total(layer_aod[indexes])
   endif else begin
       return, 0.0
   endelse
end

pro Set_Aerosol_Values, pressure_data, aerosol_data, prefix, tc_idx, summ_data, summ_srch_names
   ;; Calculate total AOD
   tot_aod = total_aod(pressure_data, aerosol_data, LAYER_AOD=layer_aod)
   total_col_name = prefix + '_TOTAL_AOD'

   where_total_column = Where(summ_srch_names eq total_col_name)
   if where_total_column[0] eq -1 then $
     Message, 'Could not find column named: ' + total_col_name

   summ_data[ where_total_column, tc_idx ] = Format_Summary_Value( tot_aod )

   ;; Calculate low
   low_index = Where(pressure_data gt 800E2, num_low)
   if low_index[0] ne -1 then begin
       low_aod = total_aod_values(layer_aod, low_index)
   endif else begin
       low_aod = 0.0
   endelse

   low_col_name = prefix + '_LOW_AOD'
   where_low_column = Where(summ_srch_names eq low_col_name)
   if where_low_column[0] eq -1 then $
     Message, 'Could not find column named: ' + low_col_name

   summ_data[ where_low_column, tc_idx ] = Format_Summary_Value( low_aod )

   ;; Calculate mid
   mid_index = Where(pressure_data gt 500E2 and pressure_data le 800E2, num_mid)
   if mid_index[0] ne -1 then begin
       mid_aod = total_aod_values(layer_aod, mid_index)
   endif else begin
       mid_aod = 0.0
   endelse

   mid_col_name = prefix + '_MID_AOD'
   where_mid_column = Where(summ_srch_names eq mid_col_name)
   if where_mid_column[0] eq -1 then $
     Message, 'Could not find column named: ' + mid_col_name

   summ_data[ where_mid_column, tc_idx ] = Format_Summary_Value( mid_aod )

   ;; Calculate high AOD
   high_index = Where(pressure_data le 500E2, num_high)
   if high_index[0] ne -1 then begin
       high_aod = total_aod_values(layer_aod, high_index)
   endif else begin
       high_aod = 0.0
   endelse

   high_col_name = prefix + '_HIGH_AOD'
   where_high_column = Where(summ_srch_names eq high_col_name)
   if where_high_column[0] eq -1 then $
     Message, 'Could not find column named: ' + high_col_name

   summ_data[ where_high_column, tc_idx ] = Format_Summary_Value( high_aod )
end

pro Testset_Summary, testcase_paths, testcase_names, summary_file, overwrite=overwrite, input_only=input_only

   ;; Used for extracting state vector elements
   control_file_base  = 'out/control1/control_file.dat'
   statevector_base   = 'statevector.dat'

   ;; For obs info/sounding info elements
   soundinginfo_base      = 'soundinginfo*.dat'
   obsinfo_base           = 'obsinfo*.nml'

   ;; Used for AOD retrieval calculation
   pressure_ret_base  = 'out/control1/pressure_levels.dat'
   pressure_in_base   = 'in/scene/*/pressure*.dat'
   aerosol_ret_sv_col  = 'AEROSOL_'

   ;; Used for getting results.dat information
   results_file_base  = 'out/results*.dat'

   sv_names_file      = 'out/aggregator/sv_names.dat'

   invalid_str = '-999.0'
  
   ;; Items to search for in files, state vector etc. Use prefixes:
   ;; FILE_ take named item(s) from gather_files_base files columns or
   ;;       headers
   ;; SND_  take named item from soundinginfo file
   ;; OBS_  take named item from obsinfo file
   ;; SV_RET_ take named item(s) from state vector retrieved column
   ;; SV_AP_ take named item(s) from state vector apriori column


   IF Keyword_Set(input_only) THEN BEGIN
       ;; Used for AOD calculation
       pressure_true_base = 'pressure*.dat'
       aerosol_true_base  = 'aerosol*.dat'
       
       gather_files_base  = [ 'psurf*.dat', 'albedo*.dat', $
                              'windspeed*.dat' $
                            ]

       summ_srch_names = [ 'Case', $

                           'CALC_TRUE_ALL_TOTAL_AOD', $
                           'CALC_TRUE_ALL_LOW_AOD', 'CALC_TRUE_ALL_MID_AOD', 'CALC_TRUE_ALL_HIGH_AOD', $

                           'CALC_TRUE_ICE_TOTAL_AOD', $
                           'CALC_TRUE_ICE_LOW_AOD', 'CALC_TRUE_ICE_MID_AOD', 'CALC_TRUE_ICE_HIGH_AOD', $

                           'CALC_TRUE_WATER_TOTAL_AOD', $
                           'CALC_TRUE_WATER_LOW_AOD', 'CALC_TRUE_WATER_MID_AOD', 'CALC_TRUE_WATER_HIGH_AOD', $

                           'CALC_TRUE_CONT_TOTAL_AOD', $
                           'CALC_TRUE_CONT_LOW_AOD', 'CALC_TRUE_CONT_MID_AOD', 'CALC_TRUE_CONT_HIGH_AOD', $

                           'CALC_TRUE_OCEANIC_TOTAL_AOD', $
                           'CALC_TRUE_OCEANIC_LOW_AOD', 'CALC_TRUE_OCEANIC_MID_AOD', 'CALC_TRUE_OCEANIC_HIGH_AOD', $

                           'SND_sounding_solar_zenith', $
                           'FILE_PSURF', $
                           'FILE_ALBEDO_1', 'FILE_ALBEDO_2', 'FILE_ALBEDO_3' ]

       summ_col_names  = [ 'Case', $

                           'CALC_TRUE_ALL_TOTAL_AOD', $
                           'CALC_TRUE_ALL_LOW_AOD', 'CALC_TRUE_ALL_MID_AOD', 'CALC_TRUE_ALL_HIGH_AOD', $

                           'CALC_TRUE_ICE_TOTAL_AOD', $
                           'CALC_TRUE_ICE_LOW_AOD', 'CALC_TRUE_ICE_MID_AOD', 'CALC_TRUE_ICE_HIGH_AOD', $

                           'CALC_TRUE_WATER_TOTAL_AOD', $
                           'CALC_TRUE_WATER_LOW_AOD', 'CALC_TRUE_WATER_MID_AOD', 'CALC_TRUE_WATER_HIGH_AOD', $

                           'CALC_TRUE_CONT_TOTAL_AOD', $
                           'CALC_TRUE_CONT_LOW_AOD', 'CALC_TRUE_CONT_MID_AOD', 'CALC_TRUE_CONT_HIGH_AOD', $

                           'CALC_TRUE_OCEANIC_TOTAL_AOD', $
                           'CALC_TRUE_OCEANIC_LOW_AOD', 'CALC_TRUE_OCEANIC_MID_AOD', 'CALC_TRUE_OCEANIC_HIGH_AOD', $


                           'SZA', $
                           'PSURF', $
                           'ALBEDO_1', 'ALBEDO_2', 'ALBEDO_3' ]

   ENDIF ELSE BEGIN
       ;; Used for AOD calculation
       pressure_true_base = 'pressure_true*.dat'
       aerosol_true_base  = 'aerosol_true*.dat'

       gather_files_base  = [ 'out/control1/x_target.dat', $
                              'x_target_true*.dat', $
                              'psurf_true*.dat', $
                              'albedo_true*.dat', $
                              'windspeed_true*.dat', $
                              'aerosol_od_true*.dat', $
                              'in/l1b/spec/OCO/*.snd', $
                              'in/l1b/spec/OCO/*.dat', $
                              'in/l1b/spec/OCO/*.0001', $
                              'in/l1b/spec/*.dat' $
                            ]

       summ_srch_names = [ 'Case', $
                           'FILE_true_x_target', $
                           'FILE_x_target', 'FILE_a_priori', 'FILE_error',  $

                           'CALC_TRUE_ALL_TOTAL_AOD', $
                           'CALC_TRUE_ALL_LOW_AOD', 'CALC_TRUE_ALL_MID_AOD', 'CALC_TRUE_ALL_HIGH_AOD', $

                           'FILE_ice_od', $

                           'CALC_TRUE_ICE_TOTAL_AOD', $
                           'CALC_TRUE_ICE_LOW_AOD', 'CALC_TRUE_ICE_MID_AOD', 'CALC_TRUE_ICE_HIGH_AOD', $

                           'FILE_water_od', $

                           'CALC_TRUE_WATER_TOTAL_AOD', $
                           'CALC_TRUE_WATER_LOW_AOD', 'CALC_TRUE_WATER_MID_AOD', 'CALC_TRUE_WATER_HIGH_AOD', $

                           'FILE_aerosol_od', $

                           'CALC_TRUE_CONT_TOTAL_AOD', $
                           'CALC_TRUE_CONT_LOW_AOD', 'CALC_TRUE_CONT_MID_AOD', 'CALC_TRUE_CONT_HIGH_AOD', $

                           'CALC_TRUE_OCEANIC_TOTAL_AOD', $
                           'CALC_TRUE_OCEANIC_LOW_AOD', 'CALC_TRUE_OCEANIC_MID_AOD', 'CALC_TRUE_OCEANIC_HIGH_AOD', $

                           'RET_ALL_TOTAL_AOD', $
                           'RET_ALL_LOW_AOD', 'RET_ALL_MID_AOD', 'RET_ALL_HIGH_AOD', $

                           'RESULTS_3', $ ; year from results.dat for FTS and OCO cases
                           'RESULTS_4', $ ; day from results.dat for FTS and OCO cases
                           'RESULTS_5', $ ; zpdtime from results.dat for FTS and OCO cases

                           'RESULTS_6', $ ; latitude from results.dat for FTS and OCO cases
                           'RESULTS_7', $ ; longitude from results.dat for FTS and OCO cases

                           'RESULTS_8', $  ; SZA 1 from results.dat for FTS and OCO cases
                           'RESULTS_18', $ ; SZA 2 from results.dat for FTS and OCO cases
                           'RESULTS_28', $ ; SZA 2 from results.dat for FTS and OCO cases

                           'FILE_sounding_solar_zenith', $ ; solar zenith from input spectrum file
                           
                           'RESULTS_10', $ ; observation zenith 1 from results.dat for FTS and OCO cases
                           'RESULTS_20', $ ; observation zenith 2 from results.dat for FTS and OCO cases
                           'RESULTS_30', $ ; observation zenith 3 from results.dat for FTS and OCO cases

                           'FILE_frame_time_stamp', $
                           'FILE_PSURF', 'SV_AP_PSURF', 'SV_RET_PSURF', $
                           'FILE_WINDSPEED', $
                           'SV_RET_WINDSPEED', $
                           'FILE_ALBEDO_1', 'FILE_ALBEDO_2', 'FILE_ALBEDO_3', $
                           'SV_RET_ALBEDO#0', 'SV_RET_ALBEDO#2', 'SV_RET_ALBEDO#4' ]

       summ_col_names  = [ 'Case', $
                           'XTARG_TRUE', $
                           'XTARG_RET', 'XTARG_RET_AP', 'XTARG_RET_VAR', $                          

                           'CALC_TRUE_ALL_TOTAL_AOD', $
                           'CALC_TRUE_ALL_LOW_AOD', 'CALC_TRUE_ALL_MID_AOD', 'CALC_TRUE_ALL_HIGH_AOD', $

                           'REPORTED_ICE_AOD', $

                           'CALC_TRUE_ICE_TOTAL_AOD', $
                           'CALC_TRUE_ICE_LOW_AOD', 'CALC_TRUE_ICE_MID_AOD', 'CALC_TRUE_ICE_HIGH_AOD', $

                           'REPORTED_WATER_AOD', $

                           'CALC_TRUE_WATER_TOTAL_AOD', $
                           'CALC_TRUE_WATER_LOW_AOD', 'CALC_TRUE_WATER_MID_AOD', 'CALC_TRUE_WATER_HIGH_AOD', $

                           'REPORTED_CONT_OCEANIC_AOD', $

                           'CALC_TRUE_CONT_TOTAL_AOD', $
                           'CALC_TRUE_CONT_LOW_AOD', 'CALC_TRUE_CONT_MID_AOD', 'CALC_TRUE_CONT_HIGH_AOD', $

                           'CALC_TRUE_OCEANIC_TOTAL_AOD', $
                           'CALC_TRUE_OCEANIC_LOW_AOD', 'CALC_TRUE_OCEANIC_MID_AOD', 'CALC_TRUE_OCEANIC_HIGH_AOD', $

                           'RET_ALL_TOTAL_AOD', $
                           'RET_ALL_LOW_AOD', 'RET_ALL_MID_AOD', 'RET_ALL_HIGH_AOD', $

                           'YEAR', $
                           'DAY', $
                           'ZPD_TIME', $

                           'LATITUDE', $
                           'LONGITUDE', $

                           'SOLAR_ZENITH_1', $
                           'SOLAR_ZENITH_2', $
                           'SOLAR_ZENITH_3', $

                           'INPUT_SZA', $

                           'OBS_ZENITH_1', $
                           'OBS_ZENITH_2', $
                           'OBS_ZENITH_3', $

                           'FRAME_TIME_STAMP', $
                           'TRUE_PSURF', 'AP_PSURF', 'RET_PSURF', $
                           'TRUE_WINDSPEED', $
                           'RET_WINDSPEED', $
                           'TRUE_ALBEDO_1', 'TRUE_ALBEDO_2', 'TRUE_ALBEDO_3', $
                           'RET_ALBEDO_1', 'RET_ALBEDO_2', 'RET_ALBEDO_3' ]
   ENDELSE

   summ_types = summ_srch_names[ uniq(summ_srch_names) ]

   if N_Elements(summ_srch_names) NE N_Elements(summ_col_names) then $
     Message, String('summ_srch_names size: ', N_Elements(summ_srch_names), $
                     ' must equal summ_col_names: ', N_Elements(summ_col_names), $
                     FORMAT='(A0, I0, A0, I0)')              

   if N_Elements(summary_file) LE 0 then $
     summary_file = 'summary.dat'

   if N_Elements(testcase_names) LE 0 then $
     testcase_names = testcase_paths

   if N_Elements(testcase_paths) NE N_Elements(testcase_names) then $
     Message, 'testcase paths and names arrays must be of the same size'

   summ_data = StrArr( N_Elements(summ_srch_names), N_Elements(testcase_paths) )
   summ_data[*] = invalid_str

   format_columns = Replicate(1, N_Elements(summ_srch_names))
   
   ;; Do not format these columns, treat as strings
   format_columns[ Where(summ_srch_names EQ 'Case') ] = 0

   skip_map = intarr(n_elements(testcase_names))
   if file_test(summary_file) and not keyword_set(overwrite) then begin
       print, 'Loading existing data from summary file: ', summary_file
       
       orig_obj  = obj_new('matrix_file', summary_file)

       if orig_obj->get_num_columns() gt 0 then begin
           orig_data = orig_obj->get_data(/NOCASTDOUBLE)
           orig_cols = orig_obj->get_all_column_labels()

           case_idx = where(orig_cols eq 'Case')
           
           for tc_idx = 0, orig_obj->get_num_rows()-1 do begin
               where_new_row = where(stregex(testcase_names, orig_data[case_idx, tc_idx], /boolean) eq 1)
               
               if where_new_row[0] ne -1 then begin
                   skip_map[where_new_row] = 1
                   
                   for col_idx = 0, orig_obj->get_num_columns()-1 do begin
                       where_new_col = Where(summ_col_names eq orig_cols[col_idx])
                       summ_data[where_new_col[0], where_new_row[0]] = orig_data[col_idx, tc_idx]
                   endfor
               endif
           
           endfor
       endif

   endif

   print, 'Parsing data for summary file: ', summary_file

   ;; Load object here so that header keywords can be added
   summ_obj = obj_new('matrix_file')
   summ_obj->Set_Column_Format_Flag, format_columns

   for tc_idx = 0, N_Elements(testcase_paths)-1 do begin
       tc_path = testcase_paths[tc_idx]
       tc_name = testcase_names[tc_idx]

       if skip_map[tc_idx] eq 1 then begin
           print, 'Existing data found, skipping testcase: ', tc_name
           continue
       endif

       print, 'Parsing testcase: ', tc_name

       ;; Add space so never interpreted as a number
       summ_data[ Where(summ_srch_names EQ 'Case'), tc_idx ] = tc_name

       for gat_idx = 0, N_Elements(gather_files_base)-1 do begin
           curr_file = File_Search(tc_path, gather_files_base[gat_idx])

           if strlen(curr_file[0]) ne 0 then begin
               curr_file_obj  = obj_new('matrix_file', curr_file[0])

               if curr_file_obj->has_valid_data() then begin
                   curr_file_data = curr_file_obj->get_data()
                   curr_file_labels = curr_file_obj->get_all_column_labels()
                   curr_file_headers = curr_file_obj->get_header_keyword_names()

                   for summ_idx = 0, N_Elements(summ_types)-1 do begin
                       col_srch_lbl = strmid( summ_types[summ_idx], $
                                              strpos(summ_types[summ_idx], 'FILE_') + 5 )

                       upper_file_labels = strupcase(curr_file_labels)
                       upper_srch_lbl = strupcase(col_srch_lbl)
                       where_curr_lbl = Where( upper_file_labels  eq upper_srch_lbl )

                       if where_curr_lbl[0] ne -1 then begin
                           summ_size = N_Elements( where(summ_srch_names eq summ_types[summ_idx]) )

                           if summ_size ne n_elements(where_curr_lbl) then $
                             Message, 'Number of matching summary types ' + summ_size + ' does not match number of matching labels: ' + n_elements(where_curr_lbl)

                           summ_data[ Where(summ_srch_names eq summ_types[summ_idx]), tc_idx ] = $
                             Format_Summary_Value( curr_file_data[ where_curr_lbl[0:summ_size-1] ] )
                       endif

                       where_curr_head = Where( curr_file_headers eq col_srch_lbl )
                       if where_curr_head[0] ne -1 then begin
                           keyword_val = curr_file_obj->get_header_keyword(col_srch_lbl)
                           key_val_parts = strsplit(keyword_val, /extract)

                           summ_data[ Where(summ_srch_names eq summ_types[summ_idx]), tc_idx ] = $
                             Format_Summary_Value( key_val_parts[0] )
                       endif

                   endfor
               endif
           endif
       endfor
       
       ;; read soundinginfo
       where_snd_present = where( stregex(summ_types, 'SND_', /BOOLEAN, /FOLD_CASE) eq 1 )
       
       if where_snd_present[0] ne -1 then begin
           soundinginfo_file = File_Search(tc_path, soundinginfo_base)
           if strlen(soundinginfo_file[0]) ne 0 then begin
               soundinginfo_vals = read_soundinginfo_file(soundinginfo_file[0], DATA_NAMES=soundinginfo_names)

               for where_idx = 0, N_Elements(where_snd_present)-1 do begin
                   summ_idx = where_snd_present[where_idx]

                   obs_srch_name = strmid( summ_types[summ_idx], $
                                           strpos(summ_types[summ_idx], 'SND_') + 4 )
                   
                   where_soundinginfo = Where( soundinginfo_names eq obs_srch_name )
                   
                   if where_soundinginfo[0] ne -1 then $
                     summ_data[ Where(summ_srch_names eq summ_types[summ_idx]), tc_idx ] = $
                       Format_Summary_Value( soundinginfo_vals[ where_soundinginfo ] )
               endfor
           endif
       endif
       
       ;; read obsinfo
       where_obs_present = where( stregex(summ_types, 'OBS_', /BOOLEAN, /FOLD_CASE) eq 1 )

       if where_obs_present[0] ne -1 then begin
           obsinfo_file = File_Search(tc_path, obsinfo_base)
           if strlen(obsinfo_file[0]) ne 0 then begin
               obsinfo_vals = read_soundinginfo_file(obsinfo_file[0], DATA_NAMES=obsinfo_names)

               for where_idx = 0, N_Elements(where_obs_present)-1 do begin
                   summ_idx = where_obs_present[where_idx]

                   obs_srch_name = 'obs_info%' + strmid( summ_types[summ_idx], $
                                                         strpos(summ_types[summ_idx], 'OBS_') + 4 )
                   
                   where_obsinfo = Where( obsinfo_names eq obs_srch_name )
                   
                   if where_obsinfo[0] ne -1 then $
                     summ_data[ Where(summ_srch_names eq summ_types[summ_idx]), tc_idx ] = $
                       Format_Summary_Value( obsinfo_vals[ where_obsinfo ] )
               endfor
           endif
       endif

       ;; read results.dat
       where_results_present = where( stregex(summ_types, 'RESULTS_', /BOOLEAN, /FOLD_CASE) eq 1 )
       
       if where_results_present[0] ne -1 then begin
           results_file = File_Search(tc_path, results_file_base)
           if strlen(results_file[0]) ne 0 then begin
               OpenR, lun_res, results_file, /GET_LUN

               ; Read header
               line_str = ''
               Readf, lun_res, line_str

               while strpos(line_str, 'Result') ge 0 and not eof(lun_res) do begin
                   Readf, lun_res, line_str
               end
              
               ; Read data
               results_str = line_str ; Last line read should be data content
               while not eof(lun_res) do begin
                   Readf, lun_res, line_str
                   if line_str ne '' then $
                     results_str = results_str + ' ' + line_str
               endwhile
               Free_Lun, lun_res
               
               ;; Make sure never read header lines as values
               if (strpos(results_str, 'Result file for') eq -1) then begin

                   results_vals = strsplit(results_str, /EXTRACT)
                   
                   for where_idx = 0, N_Elements(where_results_present)-1 do begin
                       summ_idx = where_results_present[where_idx]
                       
                       results_val_idx = fix( strmid( summ_types[summ_idx], $
                                                      strpos(summ_types[summ_idx], 'RESULTS_') + 8 ) )
                       
                       if (results_val_idx lt n_elements(results_vals)) then $
                         summ_data[ Where(summ_srch_names eq summ_types[summ_idx]), tc_idx ] = $
                           Format_Summary_Value( results_vals[ results_val_idx ] )
                   endfor
               endif 
           endif
       endif


       ;; extract required statevector columns
       statevector_file = File_Search(tc_path, statevector_base)
       if strlen(statevector_file[0]) ne 0 then begin 
           ;; Only need control file if a state vector exists
           control_file_name = File_Search(tc_path, control_file_base)
           if strlen(control_file_name[0]) eq 0 then $
             Message, control_file_base + ' not found'

           sv_file_obj = obj_new('matrix_file', statevector_file[0])
           sv_file_data = sv_file_obj->get_data()
           sv_col_lbls = sv_file_obj->get_all_column_labels()

           where_apriori = Where(sv_col_lbls eq 'A_Priori')
           if where_apriori[0] eq -1 then $
             Message, 'A_Priori column not found in ' + statevector_file

           where_sv_col = Where(sv_col_lbls eq 'Statevector')
           if where_sv_col[0] eq -1 then $
             where_sv_col = Where(sv_col_lbls eq 'State Vector')
           if where_sv_col[0] eq -1 then $
             Message, 'Statevector column not found in ' + statevector_file

           sv_col = (sv_file_data[  where_sv_col[0], * ])[*]
           sv_ap  = (sv_file_data[ where_apriori[0], * ])[*]

           ;; Apply any updates if there is such a column
           where_update_col = Where(sv_col_lbls eq 'Last Update')
           if where_update_col[0] ne -1 then $
             sv_col = sv_col + (sv_file_data[  where_update_col[0], * ])[*]

           ;; Read control file names
           names_file = File_Search(tc_path, sv_names_file)
           if strlen(names_file[0]) ne 0 then begin
               sv_names_all = Transpose((Read_Statevector_Names(names_file[0], /GROUP_AEROSOL))[*])
               sv_name_col = 0
           endif else begin
               sv_names_all = Read_Control_File_Names(control_file_name[0], SV_SIZES=sv_arr_sizes, /GROUP_AEROSOL)
               sv_name_col = where(sv_arr_sizes eq n_elements(sv_col))
           endelse
                    
           ;; Calculate  unique state vector element types
           sv_types = (sv_names_all[0,*])[ uniq(sv_names_all[0,*]) ]

           if sv_name_col[0] ne -1 then begin
               sv_names = (sv_names_all[sv_name_col[0], *])[*]
           endif else begin
               Message, 'Could not find statevector names set up matching read statevector length'
           endelse

           for summ_idx = 0, N_Elements(summ_types)-1 do begin
               if strpos(summ_types[summ_idx], 'SV_') LT 0 then $
                 continue

               for sv_type = 1, 2 do begin
                   switch sv_type of
                       1: begin
                           sv_srch_name = strmid( summ_types[summ_idx], $
                                                  strpos(summ_types[summ_idx], 'SV_RET_') + strlen('SV_RET_') )
                           sv_data = sv_col
                           break
                       end
                       2: begin
                           sv_srch_name = strmid( summ_types[summ_idx], $
                                                  strpos(summ_types[summ_idx], 'SV_AP_') + strlen('SV_AP_') )
                           sv_data = sv_ap
                           break
                       end
                   endswitch

                   sv_srch_parts = strsplit(sv_srch_name, '#', /EXTRACT)
                   sv_srch_name = sv_srch_parts[0]

                   where_sv = Where( sv_names eq sv_srch_name )

                   if where_sv[0] ne -1 and n_elements(sv_srch_parts) gt 1 then begin
                       if n_elements(where_sv) gt sv_srch_parts[1] then $
                         where_sv = where_sv[fix(sv_srch_parts[1])]
                   endif

                   if where_sv[0] ne -1 then begin
                       where_names = where(summ_srch_names eq summ_types[summ_idx])
                       summ_size = N_Elements( where_names )
                       summ_data[ where_names, tc_idx ] = $
                         Format_Summary_Value( sv_data[ where_sv[0:summ_size-1] ] )
                   endif
               endfor
           endfor
       endif

       ;; calculate total aerosol optical depth from pressure and aerosol files
       pressure_true_file = File_Search(tc_path, pressure_true_base)
       if strlen(pressure_true_file[0]) ne 0 then begin
           press_obj = obj_new('matrix_file', pressure_true_file[0])
           pressure_true_data = press_obj->get_column_data('Pressure')
       endif

       aerosol_file = File_Search(tc_path, aerosol_true_base)
       if strlen(aerosol_file[0]) ne 0 then begin
           aero_obj = obj_new('matrix_file', aerosol_file[0])
           aero_col_names = aero_obj->get_all_column_labels()
           aero_file_data = aero_obj->get_data()

           aero_locations = Where(aero_col_names ne 'Pressure' and aero_col_names ne 'pressure')

           if aero_locations[0] ne -1 then begin
               ;; Load pressure from aerosol file if not available otherwise
               where_aer_pres =  Where(aero_col_names eq 'Pressure' or aero_col_names eq 'pressure')
               if strlen(pressure_true_file[0]) eq 0 and where_aer_pres[0] ne -1 then begin
                   pressure_true_file = aerosol_file[0]
                   pressure_true_data = aero_obj->get_column_data('Pressure')
               endif

               ;; Use true psurf for bottom level
               where_true_psurf = Where(summ_col_names eq 'TRUE_PSURF')
               if where_true_psurf[0] eq -1 then $
                 Message, 'Could not find TRUE_PSURF column index'
               true_psurf = double((summ_data[ where_true_psurf, tc_idx ])[0])
               
               where_levels = where(true_psurf lt pressure_true_data)
               if where_levels[0] eq -1 then begin
                   where_levels = [ n_elements(pressure_true_data) - 1 ]
               endif

               pressure_true_data = pressure_true_data[ 0:where_levels[0] ]
               pressure_true_data[n_elements(pressure_true_data)-1] = true_psurf
               
               ;; Parse out aerosol data for specific number of levels
               aero_file_data = aero_file_data[ *, 0:where_levels[0] ]
               aerosol_data = aero_file_data[ aero_locations, * ]

               ;; Parse out aersol type data
               where_ice_col = Where(stregex(StrLowCase(aero_col_names), 'ic_?') GE 0)
               if where_ice_col[0] ne -1 then begin
                   ice_aer_data = aero_file_data[where_ice_col, * ]
               endif else begin
                   ice_aer_data = replicate(0.0, 1, n_elements(pressure_true_data))
               endelse
              
               where_water_col = Where(stregex(StrLowCase(aero_col_names), 'wc_?') GE 0)
               if where_water_col[0] ne - 1 then begin
                   water_aer_data = aero_file_data[where_water_col, * ]
               endif else begin
                   water_aer_data = replicate(0.0, 1, n_elements(pressure_true_data))
               endelse

               where_cont_col = Where(stregex(StrLowCase(aero_col_names), 'cont') GE 0)
               where_aero1_col = Where(stregex(StrLowCase(aero_col_names), 'aero1') GE 0)
               if where_cont_col[0] ne - 1 then begin
                   if where_aero1_col[0] ne -1 then begin
                       where_cont_col = [where_cont_col, where_aero1_col]
                   endif
                   cont_aer_data = aero_file_data[where_cont_col, * ]
               endif else begin
                   cont_aer_data = replicate(0.0, 1, n_elements(pressure_true_data))
               endelse

               where_oceanic_col = Where(stregex(StrLowCase(aero_col_names), 'oceanic') GE 0)
               where_aero2_col = Where(stregex(StrLowCase(aero_col_names), 'aero2') GE 0)
               if where_oceanic_col[0] ne -1 then begin
                   if where_aero2_col[0] ne -1 then begin
                       where_oceanic_col = [where_oceanic_col, where_aero2_col]
                   endif

                   oceanic_aer_data = aero_file_data[where_oceanic_col, * ]
               endif else begin
                   oceanic_aer_data = replicate(0.0, 1, n_elements(pressure_true_data))
               endelse
               
           endif
       endif

       if strlen(pressure_true_file[0]) ne 0 and strlen(aerosol_file[0]) ne 0 then begin
           ;; Calculate total AOD
           Set_Aerosol_Values, pressure_true_data, aerosol_data, 'CALC_TRUE_ALL', tc_idx, summ_data, summ_srch_names

           ;; Calculate AOD for ice clouds
           if n_elements(ice_aer_data) gt 0 then begin
               Set_Aerosol_Values, pressure_true_data, ice_aer_data, 'CALC_TRUE_ICE', tc_idx, summ_data, summ_srch_names
           endif

           ;; Calculate AOD for water clouds
           if n_elements(water_aer_data) gt 0 then begin
               Set_Aerosol_Values, pressure_true_data, water_aer_data, 'CALC_TRUE_WATER', tc_idx, summ_data, summ_srch_names
           endif

           ;; Calculate AOD for continental aerosol
           if n_elements(cont_aer_data) gt 0 then begin
               Set_Aerosol_Values, pressure_true_data, cont_aer_data, 'CALC_TRUE_CONT', tc_idx, summ_data, summ_srch_names
           endif

           ;; Calculate AOD for oceanic aerosol
           if n_elements(oceanic_aer_data) gt 0 then begin
               Set_Aerosol_Values, pressure_true_data, oceanic_aer_data, 'CALC_TRUE_OCEANIC', tc_idx, summ_data, summ_srch_names
           endif

       endif

       ;; Calculate retrieved AOD if we did find a state vector
       if n_elements(sv_names) GT 0 then begin
           where_aero_types = Where( stregex( sv_types, aerosol_ret_sv_col, /BOOLEAN ) EQ 1)

           pressure_ret_file = File_Search(tc_path, pressure_ret_base)
           if strlen(pressure_ret_file[0]) ne 0 then begin
               press_obj = obj_new('matrix_file', pressure_ret_file[0])
               pressure_ret_data = press_obj->get_column_data('Pressure')
           endif
           
           if strlen(pressure_ret_file[0]) eq 0 then begin
               pressure_in_file = File_Search(tc_path, pressure_in_base)               
               if strlen(pressure_in_file[0]) ne 0 then begin
                   press_obj = obj_new('matrix_file', pressure_in_file[0])
                   pressure_ret_data = press_obj->get_column_data('Pressure')
                   summ_obj->Set_Header_Keyword, 'Retrieval Pressure Grid', 'NOTE: Retrieval Pressure Used is Input Pressure Grid'
               endif
           endif

           if n_elements(pressure_ret_data) gt 0 and where_aero_types[0] ne -1 then begin
               aero_types = sv_types[where_aero_types]

               ;; Sizes for new aerosol matrix
               num_aero = N_Elements(aero_types)
               aero_len = N_Elements( sv_names[ Where(sv_names eq aero_types[0]) ] )
              
               aerosol_ret = DblArr( num_aero, aero_len )
               for aero_idx = 0, num_aero-1 do begin
                   where_aero_type = Where(sv_names eq aero_types[aero_idx])
                   aerosol_ret[aero_idx, *] = sv_col[ where_aero_type ]
               endfor

               ;; Check for psurf and set ground pressure to it
               where_psurf = Where( sv_names eq 'PSURF' )
               if where_psurf[0] ne -1 then begin
                   pressure_ret_data[N_Elements(pressure_ret_data)-1] = sv_col[where_psurf]
               endif

               ;; Calculate total retrieved AOD now
               Set_Aerosol_Values, pressure_ret_data, aerosol_ret, 'RET_ALL', tc_idx, summ_data, summ_srch_names
           endif
       endif
       
   endfor

   ;; Remove empty columns
   valid_col_mask = replicate( 1, n_elements(summ_col_names) )
   for col_idx = 0, n_elements(summ_col_names)-1 do begin
       where_invalid = where(summ_data[col_idx, *] eq invalid_str, invalid_count)
       if invalid_count eq N_Elements(testcase_paths) then $
         valid_col_mask[col_idx] = 0
   endfor
   where_valid_col = where(valid_col_mask eq 1)
   where_invalid_col = where(valid_col_mask eq 0)

   if where_invalid_col[0] ne -1 then begin
       empty_cols = summ_col_names[where_invalid_col]

       ;; Delete invalid data
       summ_data      = summ_data[ where_valid_col, *]
       summ_col_names = summ_col_names[ where_valid_col ]

       summ_obj->Set_Header_Keyword, 'Empty_Columns', StrJoin('"' + empty_cols + '"', " ")
   endif

   print, 'Writing summary file data to: ', summary_file

   summ_obj->set_data, summ_data
   summ_obj->set_all_column_labels, summ_col_names

   IF Keyword_Set(input_only) THEN $
     summ_obj->set_file_id, 'Input Only Testset Summary File' $
   else $ 
     summ_obj->set_file_id, 'Retrieval Testset Summary File'

   summ_obj->write, summary_file

end
