;####
;### DEPRECATED by Python routines
;###

pro Testset_Stats, testcase_paths, testcase_names, stat_file, overwrite=overwrite

   fill_value = '-999.0'

   ;; Sub-interval chi2 ranges, in wavelengths
   aband_sub_int_range = [ 0.7625, 0.7655 ]
   wco2_sub_int_range  = [ 1.602,  1.6075 ]
   sco2_sub_int_range  = [ 2.051,  2.056 ]

   out_info_base      = 'out_info*.dat'
   meas_spec_base     = 'out/rad_meas.dat'
   conv_spec_base     = 'out/rad_conv.dat'
   agg_scalar_base    = 'out/aggregator/scalar.dat'

   stat_srch_names = [ 'Case', $
                       'abs_res_rms_o2', 'abs_res_rms_weak_co2', 'abs_res_rms_strong_co2', $
                       'rel_res_rms_o2', 'rel_res_rms_weak_co2', 'rel_res_rms_strong_co2', $
                       'beg_mean_snr', 'beg_mean_snr', 'beg_mean_snr', $
                       'end_mean_snr', 'end_mean_snr', 'end_mean_snr', $
                       'meas_chi2_norm_o2', 'meas_chi2_norm_weak_co2', 'meas_chi2_norm_strong_co2', $
                       'tot_apriori_chi2', 'tot_measured_chi2', 'fc_tot_apriori_chi2', 'fc_tot_measured_chi2',$ 
                       'sub_interval_chi2_o2', 'sub_interval_chi2_weak_co2', 'sub_interval_chi2_strong_co2', $
                       'd_sigma_sq', 'iterations', 'divergences', 'outcome', 'conv_flag', 'div_flag' $
                     ]
   stat_col_names  = [ 'Case', $
                       'rms_abs_1', 'rms_abs_2', 'rms_abs_3', $
                       'rms_rel_1', 'rms_rel_2', 'rms_rel_3', $
                       'beg_mean_snr_1', 'beg_mean_snr_2', 'beg_mean_snr_3', $
                       'end_mean_snr_1', 'end_mean_snr_2', 'end_mean_snr_3', $
                       'chi2_1', 'chi2_2', 'chi2_3', $
                       'tot_apriori_chi2', 'tot_measured_chi2', 'fc_tot_apriori_chi2', 'fc_tot_measured_chi2', $
                       'sub_interval_chi2_o2', 'sub_interval_chi2_weak_co2', 'sub_interval_chi2_strong_co2', $
                       'd_sigma_sq', 'iter', 'div', 'outcome', 'conv_flag', 'div_flag' $
                     ]

   if N_Elements(stat_file) LE 0 then $
     stat_file = 'stats.dat'

   if N_Elements(testcase_names) LE 0 then $
     testcase_names = testcase_paths

   if N_Elements(testcase_paths) NE N_Elements(testcase_names) then $
     Message, 'testcase paths and names arrays must be of the same size'

   stat_data = StrArr( N_Elements(stat_col_names), N_Elements(testcase_paths) )

   stat_data[*] = fill_value

   format_columns = Replicate(1, N_Elements(stat_srch_names))
   
   ;; Do not format these columns, treat as strings
   format_columns[ Where(stat_srch_names EQ 'Case') ] = 0

   skip_map = intarr(n_elements(testcase_names))
   if file_test(stat_file) and not keyword_set(overwrite) then begin
       print, 'Loading existing data from stats file: ', stat_file
       
       orig_obj  = obj_new('matrix_file', stat_file)
       
       if orig_obj->get_num_columns() gt 0 then begin
           orig_data = orig_obj->get_data(/NOCASTDOUBLE)
           orig_cols = orig_obj->get_all_column_labels()

           case_idx = where(orig_cols eq 'Case')

           for tc_idx = 0, orig_obj->get_num_rows()-1 do begin
               where_new_row = where(stregex(testcase_names, orig_data[case_idx, tc_idx], /boolean) eq 1)

               if where_new_row[0] ne -1 then begin
                   skip_map[where_new_row] = 1

                   for col_idx = 0, orig_obj->get_num_columns()-1 do begin
                       where_new_col = Where(stat_col_names eq orig_cols[col_idx])
                       stat_data[where_new_col[0], where_new_row[0]] = orig_data[col_idx, tc_idx]
                   endfor
               endif

           endfor
       endif

   endif

   print, 'Parsing data for stats file: ', stat_file

   for tc_idx = 0, N_Elements(testcase_paths)-1 do begin
       tc_path = testcase_paths[tc_idx]
       tc_name = testcase_names[tc_idx]

       if skip_map[tc_idx] eq 1 then begin
           print, 'Existing data found, skipping testcase: ', tc_name
           continue
       endif
      
       print, 'Parsing testcase: ', tc_name

       stat_data[ Where(stat_srch_names EQ 'Case'), tc_idx ] = tc_name

       ;; read out_info to fill in stat_data values
       out_info_file = File_Search(tc_path, out_info_base)
       if strlen(out_info_file[0]) ne 0 then begin
           out_info_data = read_out_info_file(out_info_file[0], columns=out_info_cols)
           
           if out_info_data[0] ne -1 then begin
               stat_col_types = stat_srch_names[ uniq(stat_srch_names) ]

               for col_idx = 0, N_Elements(stat_col_types)-1 do begin
                   where_out_info = Where( out_info_cols eq stat_col_types[col_idx] )
                   if where_out_info[0] ne -1 then begin
                       where_stat_col = Where(stat_srch_names eq stat_col_types[col_idx])
                       ;; rescale which columns to place data into
;; since
                       ;; there are less columns available from input
;; file
                       ;; ie if for fts case does not have all 3 spec
                       if n_elements(where_stat_col) > n_elements(where_out_info) then $
                         where_stat_col = where_stat_col[0:n_elements(where_out_info)-1]

                       stat_data[where_stat_col, tc_idx ] = $
                         out_info_data[ where_out_info ]
                     endif
               endfor
           endif
       endif

       ;; read scalar to fill in stat_data values
       agg_scalar_file = File_Search(tc_path, agg_scalar_base)
       if strlen(agg_scalar_file[0]) ne 0 then begin
           agg_scalar_obj = Obj_New('matrix_file', agg_scalar_file[0])
           
           for col_idx = 0, N_Elements(stat_col_types)-1 do begin
               if agg_scalar_obj->header_has_keyword(stat_col_types[col_idx]) then begin
                   where_stat_col = Where(stat_srch_names eq stat_col_types[col_idx])

                   stat_data[where_stat_col, tc_idx ] = $
                     agg_scalar_obj->get_header_keyword(stat_col_types[col_idx])
               endif
           endfor
       endif

       ;; read measured spectra file to compute relative rms
       meas_spec_file = File_Search(tc_path, meas_spec_base)
       if strlen(meas_spec_file[0]) ne 0 then begin
           meas_spec_obj = obj_new('matrix_file', meas_spec_file)
           num_points = meas_spec_obj->get_num_rows()
           meas_rad_col = (where(strlowcase(meas_spec_obj->get_all_column_labels()) eq 'radiance'))[0]
           meas_err_col = (where(strlowcase(meas_spec_obj->get_all_column_labels()) eq 'error'))[0]
           meas_spec_data = meas_spec_obj->get_data()

           where_abs_rms = where( strmid(stat_srch_names, 0, 11) eq 'abs_res_rms' )
           where_rel_rms = where( strmid(stat_srch_names, 0, 11) eq 'rel_res_rms' )
           if (where_abs_rms[0] ne -1 and where_rel_rms[0] ne 1) then begin
               num_abs_rms = n_elements(where_abs_rms)

               norm_points = intarr(num_abs_rms)
               for pidx = 0, num_abs_rms-1 do begin
                   norm_points[pidx] = (num_points/num_abs_rms)*pidx
               endfor
               
               norm_val = double( meas_spec_data[meas_rad_col, norm_points] )

               abs_rms_str = stat_data[where_abs_rms, tc_idx ]
               where_valid_rms = where(abs_rms_str ne fill_value)
               if where_valid_rms[0] ne -1 then begin
                   abs_rms_vals = double( abs_rms_str[where_valid_rms] )
                   rel_rms_vals = abs_rms_vals / norm_val

                   where_rel_rms = where_rel_rms[where_valid_rms]               
                   stat_data[where_rel_rms, tc_idx ] = rel_rms_vals
               endif
           endif

           where_snr = where( stat_srch_names eq 'beg_mean_snr' )
           if (where_abs_rms[0] ne -1 and where_rel_rms[0] ne 1) then begin

               for bidx = 0, n_elements(where_snr)-1 do begin
                   beg_p = (num_points/num_abs_rms)*bidx
                   end_p = ((num_points/num_abs_rms)*(bidx+1))-1
                   
                   radiance_pts = meas_spec_data[meas_rad_col,beg_p:beg_p+50]
                   uncert_pts   = meas_spec_data[meas_err_col,beg_p:beg_p+50]
                   snr = radiance_pts / uncert_pts

                   stat_data[where_snr[bidx], tc_idx] = mean( snr )
               end
           endif
           
           where_snr = where( stat_srch_names eq 'end_mean_snr' )
           if (where_abs_rms[0] ne -1 and where_rel_rms[0] ne 1) then begin

               for bidx = 0, n_elements(where_snr)-1 do begin
                   beg_p = (num_points/num_abs_rms)*bidx
                   end_p = ((num_points/num_abs_rms)*(bidx+1))-1
                   
                   radiance_pts = meas_spec_data[meas_rad_col,end_p-50:end_p]
                   uncert_pts   = meas_spec_data[meas_err_col,end_p-50:end_p]
                   snr = radiance_pts / uncert_pts

                   stat_data[where_snr[bidx], tc_idx] = mean( snr )
               end
           endif
       endif
           
       ;; read convolved spectra files to compute relative rms
       conv_spec_file = File_Search(tc_path, conv_spec_base)
       if strlen(meas_spec_file[0]) ne 0 and strlen(conv_spec_file[0]) ne 0 then begin
           conv_spec_obj = obj_new('matrix_file', conv_spec_file)
           conv_spec_data = conv_spec_obj->get_data()
           conv_wvl_col = (where(strlowcase(conv_spec_obj->get_all_column_labels()) eq 'wavelength'))[0]
           conv_rad_col = (where(strlowcase(conv_spec_obj->get_all_column_labels()) eq 'radiance'))[0]
           wavelengths = conv_spec_obj->get_column_data(conv_wvl_col)

           abo2_idx = where(wavelengths ge aband_sub_int_range[0] and wavelengths le aband_sub_int_range[1])
           wco2_idx = where(wavelengths ge wco2_sub_int_range[0] and wavelengths le wco2_sub_int_range[1])
           sco2_idx = where(wavelengths ge sco2_sub_int_range[0] and wavelengths le sco2_sub_int_range[1])

           residual = meas_spec_data[meas_rad_col,*] - conv_spec_data[conv_rad_col,*]

           if abo2_idx[0] ne -1 then begin
               where_stat_col = Where(stat_srch_names eq 'sub_interval_chi2_o2') 
               tot_chi2_sum = 0.0D0
               for w_idx = 0, n_elements(abo2_idx)-1 do begin
                   r_idx = abo2_idx[w_idx]
                   tot_chi2_sum = tot_chi2_sum + (residual[r_idx]/meas_spec_data[meas_rad_col,r_idx])^2
               endfor
               stat_data[where_stat_col, tc_idx ] = tot_chi2_sum
           endif

           if wco2_idx[0] ne -1 then begin
               where_stat_col = Where(stat_srch_names eq 'sub_interval_chi2_weak_co2') 
               tot_chi2_sum = 0.0D0
               for w_idx = 0, n_elements(wco2_idx)-1 do begin
                   r_idx = wco2_idx[w_idx]
                   tot_chi2_sum = tot_chi2_sum + (residual[r_idx]/meas_spec_data[meas_rad_col,r_idx])^2
               endfor
               stat_data[where_stat_col, tc_idx ] = tot_chi2_sum
           endif

           if sco2_idx[0] ne -1 then begin
               where_stat_col = Where(stat_srch_names eq 'sub_interval_chi2_strong_co2')
               tot_chi2_sum = 0.0D0
               for w_idx = 0, n_elements(sco2_idx)-1 do begin
                   r_idx = sco2_idx[w_idx]
                   tot_chi2_sum = tot_chi2_sum + (residual[r_idx]/meas_spec_data[meas_rad_col,r_idx])^2
               endfor
               stat_data[where_stat_col, tc_idx ] = tot_chi2_sum
           endif

       endif

   endfor

   print, 'Writing stat file data to: ', stat_file

   stat_obj = obj_new('matrix_file')
   stat_obj->Set_Column_Format_Flag, format_columns
   stat_obj->set_data, stat_data
   stat_obj->set_all_column_labels, stat_col_names
   stat_obj->set_file_id, 'Testset Stats File'
   stat_obj->write, stat_file

end
