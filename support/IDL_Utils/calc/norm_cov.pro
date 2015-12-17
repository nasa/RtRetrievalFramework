function norm_cov, cov_in
   sv_len = (size(cov_in))[1]
   cov_norm = dblarr(sv_len, sv_len)
   for i = 0, sv_len-1 do begin
       for j = 0, sv_len-1 do begin
           cov_norm[i,j] = cov_in[i,j] / (sqrt(cov_in[i,i]) * sqrt(cov_in[j,j]))
       endfor
   endfor

   return, cov_norm
end 
