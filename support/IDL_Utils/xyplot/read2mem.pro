pro read2mem,nfile,fname,gmissing,header,ncol,nrow,buf
; Reads the NFILE asci data files whose names are given in FNAME(NFILE).
; Each file may be in either a "spreadsheet" or one of several Ames formats.
; Entire contents of data files are stored in a contiguous 1-D array BUF.
;
;  INPUTS:
;          nfile                    Number of input files
;          fname(nfile)             Names (full path) of input files
;          gmissing                 Value to be used to represent missing data
;
; OUTPUTS:
;          header(nfile,1000)       2-D string array of column titles
;          ncol(nfile)              Number of columns in each file
;          nrow(nfile)              Number of rows in each file
;          buf(nele)                Array of data values
;
;   NOTES:
; 1)All data files are expected to contain NLHEAD,NFMT on their first line,
;   where NLHEAD is the number of header lines,
;   and NFMT is the data format code.
;   If (NFMT le 1000) the file is assumed to be in a spreadsheet format
;   with NFMT being the number of columns of data.
;   Otherwise it is assumed to be in an Ames format with NFMT representing
;   the Ames format code. This means that this subroutine cannot handle a
;   spreadsheet file having more than 1000 data columns.
; 2)Missing data values are converted from their local values, which differ from
;   file to file and parameter to parameter, to GMISSING, the global missing value.
;
   ccc = string(' ')
   iv1 = string(' ')
   iv2 = string(' ')
   kele=long(0)
   for ifile=0,nfile-1 do begin
   nfmt=long(0)
   nx=long(0)
   ndcol=long(0)
   nauxv=long(0)
   nauxc=long(0)
   openr, unit, fname(ifile), /get_lun

   readf, unit, ccc
   IF strpos(StrLowCase(ccc), 'begin ') GE 0 THEN BEGIN
       matrix_obj = Obj_New('MATRIX_FILE', fname(ifile))

       col_labels = matrix_obj->Get_All_Keyword_Items('Labels')
       header[ifile, 0:n_elements(col_labels)-1] = col_labels
       
       ncol[ifile] = matrix_obj->Get_Num_Columns()
       nrow[ifile] = matrix_obj->Get_Num_Rows()

       mat_data = (matrix_obj->Get_Data())[*]

       buf[kele:kele+n_elements(mat_data)-1] = mat_data

       kele = kele + n_elements(mat_data)

       close, unit
       free_lun, unit
       continue
   ENDIF ELSE BEGIN
       reads, ccc, nlhead, nfmt
   ENDELSE

   case 1 of

   (nfmt eq 1001): begin     ; It is an Ames format
      for k=1,7 do  readf, unit, ccc
      readf, unit, iv1
      readf, unit, ndcol
      zscal = fltarr(ndcol)
      zmiss = fltarr(ndcol)
      ztitl = strarr(ndcol)
      readf, unit, zscal
      readf, unit, zmiss
      readf, unit, ztitl
      zscal=[1.0,zscal]
      zmiss=[-999,zmiss]
      ztitl=[iv1,ztitl]
      ndcol=ndcol+1  ; the IV is grouped with the data variables
      end
   (nfmt eq 1010): begin     ; It is an Ames format
      for k=1,7 do  readf, unit, ccc
      readf, unit, iv1
      readf, unit, ndcol
      zscal = fltarr(ndcol)
      zmiss = fltarr(ndcol)
      ztitl = strarr(ndcol)
      readf, unit, zscal
      readf, unit, zmiss
      readf, unit, ztitl
;
      readf, unit, nauxv
      ascal=fltarr(nauxv)
      amiss=fltarr(nauxv)
      atitl=strarr(nauxv)
      readf, unit, ascal
      readf, unit, amiss
      readf, unit, atitl
;
      ascal=[1.0,ascal]
      amiss=[-999,amiss]
      atitl=[iv1,atitl]
      nauxv=nauxv+1  ; the IV is grouped with the auxiliary variables
      end
   (nfmt eq 2010):  begin     ; It is in an Ames format
      for k=1,7 do  readf, unit, ccc
      readf, unit, nx
      readf, unit, nxdef
      xx = strarr(nxdef)
      readf, unit, ccc
      substr,ccc,xx,nxdef,jj
      readf, unit, iv1
      readf, unit, iv2
      readf, unit, kdcol
      dscal = fltarr(kdcol)
      dmiss = fltarr(kdcol)
      dtitl = strarr(kdcol)
      readf, unit, dscal
      readf, unit, dmiss
      readf, unit, dtitl
      ndcol=kdcol*nx
      zscal=fltarr(ndcol)
      zmiss=fltarr(ndcol)
      ztitl=strarr(ndcol)
      for i=0,kdcol-1 do begin
        zscal(i*nx:i*nx+nx-1)=dscal(i)
        zmiss(i*nx:i*nx+nx-1)=dmiss(i)
        ztitl(i*nx:i*nx+nx-1)=dtitl(i)+'_'+xx
      endfor
      ascal=[1.0]
      amiss=[-999]
      atitl=[iv2]
      readf, unit, nauxv
      if nauxv gt 0. then begin
        ascal=fltarr(nauxv)
        amiss=fltarr(nauxv)
        atitl=strarr(nauxv)
        readf, unit, ascal
        readf, unit, amiss
        readf, unit, atitl
        ascal=[1.0,ascal]
        amiss=[-999,amiss]
        atitl=[iv2,atitl]
      endif
      nauxv=nauxv+1  ; IV2 is grouped with the auxiliary variables
      end
   (nfmt eq 2110):  begin     ; It is in  an Ames format
      for k=1,7 do  readf, unit, ccc
      readf, unit, iv1
      readf, unit, iv2
      readf, unit, ndcol
      zscal = fltarr(ndcol)
      zmiss = fltarr(ndcol)
      ztitl = strarr(ndcol)
      readf, unit, zscal
      readf, unit, zmiss
      readf, unit, ztitl
      zscal=[1.0,zscal]
      zmiss=[-999,zmiss]
      ztitl=[iv1,ztitl]
      readf, unit, nauxv
      ascal=fltarr(nauxv)
      amiss=fltarr(nauxv)
      atitl=strarr(nauxv)
      readf, unit, ascal
      readf, unit, amiss
      readf, unit, atitl
      ascal=[1.0,ascal]
      amiss=[-999,amiss]
      atitl=[iv2,atitl]
      ndcol=ndcol+1  ; IV1 is grouped with the data variables
      nauxv=nauxv+1  ; IV2 is grouped with the auxiliary variables
      end
   (nfmt eq 2160): begin     ; It is an Ames format
      for k=1,8 do  readf, unit, ccc
      readf, unit, iv1
      readf, unit, iv2
      readf, unit, ndcol
      zscal = fltarr(ndcol)
      zmiss = fltarr(ndcol)
      ztitl = strarr(ndcol)
      readf, unit, zscal
      readf, unit, zmiss
      readf, unit, ztitl
;
      readf, unit, nauxv
      readf, unit, nauxc
      nauxv=nauxv-nauxc
      ascal=fltarr(nauxv)
      amiss=fltarr(nauxv)
      atitl=strarr(nauxv)
      readf, unit, ascal
      readf, unit, amiss
      readf, unit, ccc                     ;  skip LENA(a)
      for j=1,nauxc do  readf, unit, ccc   ;  skip AMISS(a)
      readf, unit, atitl
      for j=1,nauxc do  readf, unit, ccc   ; skip character variable titles
;
      zscal=[1.0,zscal]
      zmiss=[-999,zmiss]
      ztitl=[iv1,ztitl]
      ndcol=ndcol+1  ;  IV1 is grouped with the data variables
;                       IV2 is a character variable and is skipped
      end
   (nfmt le 1000): begin    ;  File is in a spreadsheet format
      ndcol=nfmt
      ztitl = strarr(ndcol)
      zscal = fltarr(ndcol)
      zmiss = fltarr(ndcol)
      jj=0
      zscal(*)=1.0       ; default value
      zmiss(*)=-999.     ; default value
      for k=2,nlhead do begin
         readf, unit, ccc
         lcolon=strpos(ccc,":")+1
      if(strpos(ccc,"MISSING:") ge 0)then reads,strmid(ccc,lcolon,9999),zmiss
      if(strpos(ccc,"SCALING:") ge 0)then reads,strmid(ccc,lcolon,9999),zscal
      if(strpos(ccc,"LABELS:") ge 0)then substr,strmid(ccc,lcolon,9999),ztitl,ndcol,jj
      endfor
      if (jj eq 0) then  substr,ccc,ztitl,ndcol,jj
      if jj ne ndcol then begin
         print, fname(ifile)+': mismatched number of columns & ztitl'
         print, 'Stated # columns =',ndcol,'  Actual # =',jj
         for j=0,ndcol-1  do begin
           print,j,ztitl(j)
         endfor
      endif
      end
   else: print,'unrecognized AMES format'
   endcase
;
;
   if (nfmt gt 1000) then begin  ;  Read Ames format comment lines
      readf, unit, nscomm  ;  Skip comment lines
      for k=1,nscomm do readf, unit, ccc
      readf, unit, nncomm  ;  Skip comment lines
      for k=1,nncomm do readf, unit, ccc
   endif
   ntcol=ndcol+nauxv
   ztitl=strtrim(ztitl,2)
   print,nauxv,nauxc,ndcol,ntcol
   for icol=0,ndcol-1 do begin  ; remove any trailing junk from column labels
      lens=strpos(ztitl(icol),' ')
      if(lens gt 0) then  ztitl(icol)=strmid(ztitl(icol),0,lens)
   endfor
   for icol=0,nauxv-1 do begin  ; remove any trailing junk from column labels
      lens=strpos(atitl(icol),' ')
      if(lens gt 0) then  atitl(icol)=strmid(atitl(icol),0,lens)
   endfor
;
;
; Read through data files to get VAL and ERR arrays.
; Note that entire files are read to determine NROW, even if PNAME was not found
      if (nauxv gt 0) then  aval = fltarr(nauxv)
      zval = fltarr(ndcol)
      nrec=long(1)
      krow=long(0)
      while not EOF(unit)  do begin
         on_ioerror, closeunit
         if nauxv gt 0 then begin
            if (nfmt eq 2160) then readf, unit, ccc ; IV2=X(m,2) char string
            readf, unit, aval  ; read the auxiliary data values
;            if (nfmt eq 2160) then nrec=aval(0)
            if (nfmt eq 2160) then nrec=aval(nauxv-nauxc)  ; kluge for ILAS
            if (nfmt eq 2110) then nrec=aval(1)
            for j=1,nauxc  do readf, unit, ccc   ; skip character values
            bad=where(aval eq amiss, nmiss)
            aval=aval*ascal
            if (nmiss gt 0) then aval(bad) = gmissing
         endif
         for irec=1,nrec do begin
             readf,unit,zval
             bad=where(zval eq zmiss, nmiss)
             zval=zval*zscal
             if (nmiss gt 0) then zval(bad) = gmissing
             if nauxv gt 0 then buf(kele)=aval
             buf(kele+nauxv)=zval
             kele=kele+ntcol
             krow=krow+1
         endfor  ;  irec=1,nrec
;         krow = krow + nrec  ; nrec is frequently wrong
      endwhile
      closeunit:
      close, unit
      free_lun, unit
      if nauxv gt 0 then header(ifile,0:nauxv-1)=atitl
      header(ifile,nauxv:ntcol-1)=ztitl
      ncol(ifile) = ntcol
      nrow(ifile)=krow
   endfor  ; ifile=0,nfile-1
end
