PRO hist_l2_run_stats, list_files, summ_files, brdf_map_file, set_dir, plot_file

; DEPENDENCIES:
; ** Locations give in parantheses. Root directory is /home/odell/idlprogs/
;
;   barplot2  (libs/colib/)
;       barplot (libs/colib/)
;   legend (libs/printlib/)


NUM_BINS = 9
BIN_WIDTH = 5

oldp=!p
oldx=!x
oldy=!y

tvlct, ct_r, ct_g, ct_b, /get
loadct, 39

!p.multi = [0,2,1]
!p.charsize = 1
!y.omargin = [1,0]

set_plot, 'PS'
device, /color, xsize=11,ysize=5.5, /inches, $
        xoff=1.5, yoff=11, bits = 8, /landscape, $
        filename=plot_file, set_font = 'Times-Bold',/bold	

;; Load list files
list_objs  = objarr(n_elements(list_files))
list_names = strarr(n_elements(list_files))
for list_idx = 0, N_Elements(list_files)-1 do begin
    list_base = file_basename(list_files[list_idx])
    beg_pos = 0
    end_pos = strpos(list_base, '.', /reverse_search)
    list_names[list_idx] = strmid(list_base, beg_pos, end_pos-beg_pos)

    list_objs[list_idx] = read_column_file(list_files[list_idx], CHUNKSIZE=2000L)
endfor

;; Load stat files
for summ_idx = 0, N_Elements(summ_files)-1 do begin
    summ_obj = obj_new('matrix_file', summ_files[summ_idx])
    summ_arr = summ_obj->get_struct_array()

    if n_elements(summ_data) eq 0 then begin
        summ_data = summ_arr
    endif else begin
        if n_elements(tag_names(summ_data)) gt n_elements(tag_names(summ_arr)) then begin
            summ_data_new = replicate(summ_data[0], n_elements(summ_data) + n_elements(summ_arr))
        endif else begin
            summ_data_new = replicate(summ_arr[0], n_elements(summ_data) + n_elements(summ_arr))
        endelse

        for dst_idx = 0, n_elements(summ_data)-1 do begin
            dst_struct = summ_data_new[dst_idx]
            struct_assign, summ_data[dst_idx], dst_struct
            summ_data_new[dst_idx] = dst_struct
        endfor
        src_idx = 0
        for dst_idx = n_elements(summ_data), n_elements(summ_data_new)-1 do begin
            dst_struct = summ_data_new[dst_idx]
            struct_assign, summ_arr[src_idx], dst_struct
            summ_data_new[dst_idx] = dst_struct
            src_idx = src_idx + 1
        endfor
        summ_data = summ_data_new
        summ_data_new = 0L
    endelse
endfor

;; Convert to long 64 bit integers
summ_ids = ulong64(summ_data.case_column)

;; Load brdf map file
brdf_obj = read_column_file(brdf_map_file, CHUNKSIZE=2000L)
brdf_obj->set_column_label, 0, 'case_column'
brdf_obj->set_column_label, 1, 'brdf_type'
brdf_vals = brdf_obj->get_struct_array()
brdf_types = different(brdf_vals.brdf_type)
brdf_ids = ulong64(brdf_vals.case_column)

id = intarr(n_elements(list_objs))
for brdf_idx = 0,n_elements(brdf_types)-1 do begin
    sza_bins = findgen(NUM_BINS)*10. + BIN_WIDTH
    y_bin = fltarr(NUM_BINS, n_elements(list_objs))
    number_cases = 0

    for obj_idx = 0, n_elements(list_objs)-1 do begin
        list_ids_all = ulong64((list_objs[obj_idx])->get_data(/NOCASTDOUBLE))

        is_brdf_type = replicate(0, n_elements(list_ids_all))
        for list_idx = 0, n_elements(list_ids_all)-1 do begin
            where_brdf_type = where(list_ids_all[list_idx] eq brdf_ids)

            if where_brdf_type[0] eq -1 then $
              Message, 'Could not find brdf type for item ' + strcompress(string(list_ids_all[list_idx]))
            if brdf_types[brdf_idx] eq (brdf_vals[where_brdf_type[0]]).brdf_type then $
              is_brdf_type[list_idx] = 1
        endfor

        brdf_list_indexes = where(is_brdf_type eq 1)
        if brdf_list_indexes[0] eq -1 then $
          continue

        list_ids_brdf = list_ids_all[brdf_list_indexes]

        where_data = replicate(-1L, n_elements(list_ids_brdf))
        for list_idx = 0, n_elements(list_ids_brdf)-1 do begin
            where_name = where(list_ids_brdf[list_idx] eq summ_ids)

            if where_name[0] eq -1 then $
              Message, 'Could not find data for list item ' + strcompress(string(list_ids_brdf[list_idx]))
            if n_elements(where_name) gt 1 then $
              Message, 'Too many data rows ' + strcompress(string(n_elements(where_name))) + ' match list name ' + strcompress(string(list_ids_brdf[list_idx]))

            where_data[list_idx] = where_name[0]
            number_cases = number_cases + 1
        endfor

        for bin_idx = 0, NUM_BINS-1 do begin
            where_sza_range = where(abs(summ_data[where_data].input_sza - sza_bins[bin_idx]) LE BIN_WIDTH, bin_count )
            y_bin[bin_idx,obj_idx] = bin_count
        endfor
    endfor

    if number_cases gt 0 then begin

        if size(y_bin, /n_dimensions) eq 1 then $
          ytot = total(y_bin) $
        else $
          ytot = total(y_bin,2)
    
        yr = [0,max(ytot)*1.1]
        tit = brdf_types[brdf_idx] + ', total number cases = ' + strcompress(string(number_cases))
        
        barplot2, sza_bins, y_bin, xr=[0,90], wid = .03, col = colors, yr=yr, $
                  tit=tit, xtit = 'Solar Zenith Angle [deg]', ytit = '# cases'
    endif

    if (brdf_idx eq 1 or n_elements(brdf_types) eq 1) then legend, list_names, col = colors, th=2+id+2, /top, /right, lines=0+id
endfor

xyouts, /norm, 0.5, 0.01, align=0.5, set_dir, charsize=1.5 -0.5

device, /close_file

!p = oldp
!x = oldx
!y = oldy
tvlct, ct_r, ct_g, ct_b

END
    
    



