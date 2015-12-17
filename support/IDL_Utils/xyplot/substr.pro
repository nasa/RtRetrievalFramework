     pro substr,inputstring,outputarray,mss,nss
;  Converts a space/tab/comma-delimited character string (INPUTSTRING),
;  into an array of its NSS sub-string components (OUTPUTARRAY).
;
;  Inputs:
;       INPUTSTRING    The character string to be parsed.
;       MSS            Declared dimension of OUTPUTARRAY in calling prog.
;
;   Outputs:
;       NSS           Number of sub-strings that were found in INPUTSTRING.
;       OUTPUTARRAY   Array containing the NSS sub-strings that were found.
;
;  Special Notes:
;  1) If the actual number of sub-strings exceeds MSS, only the first
;     MSS sub-strings will be copied to OUTPUTARRAY, avoiding the
;     possibility of array-bound violations.
;  3) Currently recognized delimiters include
;                       nul            (ASCII character # 0)
;                       horizontal tab (ASCII character # 9)
;                       space          (ASCII character # 32)
;                       comma          (ASCII character # 44)
;  4) Calls external functions FNBC.PRO (First Non-Blank Character),
;     and FBC.PRO (First Blank Character)
;
;  18-May-98  GCT
;
     lenin = strlen(strtrim(inputstring))
     nss=0
     iend = -1 ; ending index of substring +1
     fnbc,inputstring,ibeg
     while ibeg gt iend  do begin
        fbc,strmid(inputstring,ibeg+1,lenin-ibeg),nc
        iend=ibeg+nc+1
        if nss lt mss then outputarray(nss)=strmid(inputstring,ibeg,iend-ibeg)
        fnbc,strmid(inputstring,iend+1,lenin-iend),nc
        ibeg=iend+nc+1
        nss=nss+1
      endwhile
      return
      end
