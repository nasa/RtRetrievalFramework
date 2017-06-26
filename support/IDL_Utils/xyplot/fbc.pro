      pro fbc,label,nc
; Returns position of First Blank Character in string LABEL
; Recognized "blank" characters include:
;       Null              ASCII character # 00
;       Horizontal Tab    ASCII character # 09
;       Space             ASCII character # 32
;       Comma             ASCII character # 44
;
      for nc=0,strlen(label)-1 do begin
         cc=strmid(label,nc,1)
         if (    cc eq string(00B)   $
             or  cc eq string(09B)   $
             or  cc eq string(32B)   $
             or  cc eq string(44B) ) then return  ; Successful Return
      endfor
      return       ;  Abnormal return: No "blanks" found; NC = strlen(label)
      end
