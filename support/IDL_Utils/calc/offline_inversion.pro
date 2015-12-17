function find_sing_idx, mat_w, tol
sing_idx = -1L
mach_cond  = (machar(/double)).xmin
s_sum_corr = 0.0D0 
s_sum_m1   = 0.0D0 
tol = 99
print, 'sing_idx, s_sum_m1, s_sum_corr, mach_cond, eps' 

min_val = (machar(/double)).xmax
min_idx = -1L
for sing_idx = 0, n_elements(mat_w)-1 do begin
    if (mat_w[sing_idx] ne 0.0E0) then begin
        s_sum_m1   = s_sum_m1 + 1.0/mat_w[sing_idx]
    endif

    ;; Sum singular values explicitly so we
    ;; can check for a 0 denominator for inverse 
    s_sum_corr = 0.0D0
    for rest_idx = sing_idx, n_elements(mat_w)-1 do begin
        s_sum_corr = s_sum_corr + mat_w[rest_idx] 
    endfor
    
    mach_cond = s_sum_corr * s_sum_m1 ; * (machar(/double)).eps
    print, sing_idx, s_sum_m1, s_sum_corr, mach_cond

    if (mach_cond lt min_val) then begin
        min_val = mach_cond
        min_idx = sing_idx
    endif
endfor
print, 'min_val = ', min_val
print, 'min_idx = ', min_idx

sing_idx = sing_idx - 1

return, sing_idx
end


base_dir = '/home/mcduffie/oco_l2/testing/testcases/ecmwf_pd/extended_mixed_tests/sublayer_od_p_20/retrieval_set_01_lev_20/ECMWF_T00877/sav.f90-5628M.dflt/'
;base_dir = '/home/mcduffie/oco_l2/L2_Tests/test_nadir/pf_Jul1_TROP_01/Retrieval/run.f90-5867M'

;base_dir = '/home/mcduffie/oco_l2/testing/testcases/ecmwf_pd/extended_mixed_tests/sublayer_od_p_20/retrieval_set_01_lev_20/ECMWF_T00877/run.f90-5628M/'
;base_dir = '/home/mcduffie/oco_l2/testing/testcases/ecmwf_pd/extended_mixed_tests/sublayer_od_p_20/retrieval_set_01_lev_20/ECMWF_T00877/run.f90-5867M/'

iter_str = '.iter02'
;gamma_lm = 0.0
gamma_lm = 10.0
;gamma_lm = 5.0
;gamma_lm = 2.5

statevector_f = base_dir + '/out/statevector.dat' + iter_str
cov_ap_f = base_dir + '/out/cov_ap.dat' + iter_str
cov_ap_m1_f = base_dir + '/out/cov_ap_m1.dat' + iter_str
rad_model_f = base_dir + '/out/rad_conv.dat' + iter_str
rad_meas_f = base_dir + '/out/rad_meas.dat' + iter_str
pd_f = base_dir + '/out/pd.dat' + iter_str

if (n_elements(no_reread) EQ 0) then begin
    cov_ap_o = obj_new('matrix_file', cov_ap_f)
    cov_ap_m1_o = obj_new('matrix_file', cov_ap_m1_f)
    statevector_o = obj_new('matrix_file', statevector_f)
    rad_model_o = obj_new('matrix_file', rad_model_f)
    rad_meas_o = obj_new('matrix_file', rad_meas_f)
    pd_o = obj_new('matrix_file', pd_f)

    sv_ap = (statevector_o->get_data())[1,*]
    statevector = (statevector_o->get_data())[2,*]
    cov_ap = cov_ap_o->get_data()
    cov_ap_m1_chk = cov_ap_m1_o->get_data()
    rad_model = (rad_model_o->get_data())[2, *]
    rad_meas = (rad_meas_o->get_data())[2:3, *]
    pd = pd_o->get_data()
endif


print, 'Check that gamma_lm is really', gamma_lm, ' for ', iter_str

sv_len = n_elements(statevector)
tot_pixel = (size(rad_meas))[2]

sigma_ap = dblarr(sv_len)
for sv_idx = 0, sv_len-1 do begin
    sigma_ap[sv_idx] = sqrt(cov_ap[sv_idx, sv_idx])
endfor

statevector = statevector / sigma_ap
sv_ap = sv_ap / sigma_ap

pd_zeta = dblarr(sv_len, tot_pixel)
for i = 0, sv_len-1 do begin
    for j = 0, tot_pixel-1 do begin
        pd_zeta[i,j] = pd[i,j] * sigma_ap(i)
    endfor
endfor
pd = pd_zeta

cov_rad_m1 = 1.D0 / rad_meas[1,*]^2

cov_ap_norm = dblarr(sv_len, sv_len)
for i = 0, sv_len-1 do begin
    for j = 0, sv_len-1 do begin
        cov_ap_norm[i,j] = cov_ap[i,j] / (sqrt(cov_ap[i,i]) * sqrt(cov_ap[j,j]))
    endfor
endfor

; (n,1)=(n,n)*(n,1)
delta_ap = (statevector - sv_ap)

; compute cov_ap_m1 and cov_ap_m1_sqrt
la_svd, cov_ap_norm, mat_w1_m1, mat_u_m1, mat_v_m1

; Compute cov_ap_m1 from SVD 
sing_idx = find_sing_idx(mat_w1_m1, 1e-10)
;check_idx = sing_idx
mat_w2_m1 = dblarr(n_elements(mat_w1_m1))
wl1_last = 0.0D0
wl1_old = 0.0D0
for check_idx = sing_idx, sv_len-1 do begin
    mat_w2_m1[*] = 0.0D0
    mat_w2_m1[0:check_idx] = (1.0D0 / mat_w1_m1)[0:check_idx]
    cov_ap_m1 = mat_v_m1 ## diag_matrix(mat_w2_m1) ## transpose(mat_u_m1)

    w1 = cov_ap_m1 ## delta_ap
    wl1_old = total(transpose(delta_ap) * w1)
    ;; (x-xa)T(Sa-1)(x-xa)
    ;; Compute dot_product


    if check_idx ne sing_idx then begin
        cf_ratio = abs(wl1_old - wl1_last)/wl1_old
        print, check_idx, wl1_old, wl1_last, cf_ratio
        if wl1_old lt 0.0 or cf_ratio gt 0.2 then begin
            cov_ap_m1 = cov_ap_m1_last
            wl1_old = wl1_last
            break
        endif
    endif else begin
        print, check_idx, wl1_old
    endelse
    cov_ap_m1_last = cov_ap_m1
    wl1_last = wl1_old
endfor

la_svd, cov_ap_norm, mat_w1_m1_sqrt, mat_u_m1_sqrt, mat_v_m1_sqrt
mat_w2_m1_sqrt = dblarr(n_elements(mat_w1_m1_sqrt))
mat_w2_m1_sqrt = (1.0D0 / sqrt(mat_w1_m1_sqrt))
cov_ap_m1_sqrt = mat_v_m1_sqrt ## diag_matrix(mat_w2_m1_sqrt) ## transpose(mat_u_m1_sqrt)

; (n,m) * (m,m)
wm1 =  transpose(pd) ## diag_matrix(cov_rad_m1)

; (n,n)=scalar*(n,n)
wm3 = (1 + gamma_lm) * cov_ap_m1

; (n,n)=(n,m)*(m,n)
ktsek = wm1 ## pd

; (n,n)
wm2 = wm3 + ktsek

;;;;
;delta_ap_copy = delta_ap
;print, 'gamma2_old', 'gamma2_new', 'tt1', 'wl1_old', 'wl1_new', format='(5(A16))'
;openw, lun, '~/cf_single_diag.txt', /get_lun
;printf, lun, 'idx', 'sv-sv_ap', 'gamma2_old', 'gamma2_new', 'tt1', 'wl1_old', 'wl1_new', format='(A4,6(A16))'
;for d_idx = 0, n_elements(delta_ap)-1 do begin
;delta_ap_new = dblarr( n_elements(delta_ap) )
;delta_ap_new[d_idx] = delta_ap_copy[d_idx]
;delta_ap = delta_ap_new
;;;;

wl = (cov_ap_m1_sqrt ## delta_ap)
wl1_new = total(wl*wl)

;(y-F(x))T(Se-1)(y-F(X))
;     tt=matmul(cov_rad_m1,(rad_meas(:,1) - rad_model))
;     (m,1)=(m,m)*(m,1)
tt = dblarr(tot_pixel)
for i = 0, tot_pixel-1 do begin
    tt[i] = cov_rad_m1[i] * (rad_meas[0,i] - rad_model[i])
endfor
tt1 = total( transpose(rad_meas(0,*) - rad_model) * tt )

gamma2_old=tt1+wl1_old
gamma2_new=tt1+wl1_new

;;;;
;printf, lun, d_idx, delta_ap[d_idx], gamma2_old, gamma2_new, tt1, wl1_old, wl1_new, format='(I4,6(E16.8))'
;endfor
;free_lun, lun
;;;;

w2 = wm1 ## (rad_meas[0,*] - rad_model) - w1

; now solve wm2 * delta_sv = w2 for delta_sv
delta_sv_old = la_least_squares(wm2, w2, /double, method=0)
;delta_sv_old = la_least_squares(wm2, w2, /double, method=2, rank=out_rank, rcondition=1e-10)

la_svd, wm2, mat_w1_wm2, mat_u_wm2, mat_v_wm2

wm2_m1 = mat_v_wm2 ## diag_matrix(1.0/mat_w1_wm2) ## transpose(mat_u_wm2)
delta_sv_est = wm2_m1 ## w2

incp_sum = 0.0
min_idx = -1
min_val = (machar(/double)).xmax
check_vec = transpose(sigma_ap)
;for mse_idx = 0, sv_len-1 do begin
for mse_idx = sv_len+1, sv_len-1 do begin
    if abs(mat_w1_wm2[mse_idx]) gt 0.0 then $
      incp_sum = incp_sum + (1.0 / mat_w1_wm2[mse_idx]^2)

    decr_sum = 0.0
    for sum_idx = mse_idx+1, sv_len-1 do begin
        decr_sum = decr_sum + (mat_v_wm2[*, sum_idx]##check_vec)^2
;        decr_sum = decr_sum + ((transpose(mat_v_wm2))[*,sum_idx]##delta_sv_est)^2
    endfor

    test_sum = sigma_ap[mse_idx]^2 * incp_sum + decr_sum
;    test_sum = incp_sum + decr_sum

    print, mse_idx, sigma_ap[mse_idx]^2 * incp_sum, decr_sum, test_sum, format='(I4,4(E16.8))'
    if test_sum lt min_val then begin
        min_val = test_sum
        min_idx = mse_idx
    endif 
endfor
;print, 'min_val = ', min_val
;print, 'min_idx = ', min_idx
;check_idx = min_idx
sing_idx = find_sing_idx(mat_w1_wm2, 1e-6)
check_idx = sing_idx

mat_w2_wm2 = dblarr(n_elements(mat_w1_wm2))
;for check_idx = sing_idx, sv_len-1 do begin
    mat_w2_wm2[*] = 0.0
    mat_w2_wm2[0:check_idx] = (1.0/mat_w1_wm2)[0:check_idx]
    wm2_m1 = mat_v_wm2 ## diag_matrix(mat_w2_wm2) ## transpose(mat_u_wm2)
    delta_sv_new = wm2_m1 ## w2
   
    d_sigma_sq = total(delta_sv_new * w2)

;    if check_idx ne sing_idx then begin
;        ds_ratio = abs(d_sigma_sq_last - d_sigma_sq)/d_sigma_sq
;        print, check_idx, d_sigma_sq, d_sigma_sq_last, ds_ratio
;        if d_sigma_sq lt 0.0 or ds_ratio gt 0.2 then begin
;            d_sigma_sq = d_sigma_sq_last
;            delta_sv_new = delta_sv_new
;            break
;        endif
;    endif else begin
;        print, check_idx, d_sigma_sq
;    endelse
;    d_sigma_sq_last = d_sigma_sq
;    delta_sv_last = delta_sv_new
;endfor

la_svd, wm2, mat_w, mat_u, mat_v
mat_w = sqrt(mat_w)
wm2_sqrt = mat_v ## diag_matrix(mat_w) ## transpose(mat_u)

wml2 = wm2_sqrt ## delta_sv_new
;d_sigma_sq_new = total(transpose(delta_sv_new) * wml2)
d_sigma_sq_new = total( wml2*wml2 )

print, 'gamma2_old', 'gamma2_new', 'tt1', 'wl1_old', 'wl1_new', format='(5(A16))'
print, gamma2_old, gamma2_new, tt1, wl1_old, wl1_new

print, 'd_sigma_sq, d_sigma_sq_new'
print, d_sigma_sq, d_sigma_sq_new

end
