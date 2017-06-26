      pro fnbc,label,nc
; Returns position of First Non-Blank Character in string LABEL
; Recognized "blank" characters include:
;       Null              ASCII character # 00
;       Horizontal Tab    ASCII character # 09
;       Space             ASCII character # 32
;       Comma             ASCII character # 44
; NC=-1 indicates that the entire string was blank.
;
      for nc=0,strlen(label)-1 do begin
         cc=strmid(label,nc,1)
         if (    cc ne string(00B)   $
            and  cc ne string(09B)   $
            and  cc ne string(32B)   $
            and  cc ne string(44B) ) then return
      endfor
      nc=-1
      return       ;  Abnormal return: No characters found
      end
